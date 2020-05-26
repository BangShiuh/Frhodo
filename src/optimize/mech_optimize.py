# This file is part of Frhodo. Copyright © 2020, UChicago Argonne, LLC
# and licensed under BSD-3-Clause. See License.txt in the top-level 
# directory for license and copyright information.

from qtpy.QtCore import QThreadPool
import numpy as np
import multiprocessing as mp
from copy import deepcopy

from optimize.optimize_worker import Worker
from optimize.fit_fcn import initialize_parallel_worker, update_mech_coef_opt

class Multithread_Optimize:
    def __init__(self, parent):
        self.parent = parent
        
        # Initialize Threads
        parent.optimize_running = False
        # parent.threadpool = QThreadPool()
        # parent.threadpool.setMaxThreadCount(2) # Sets thread count to 1 (1 for gui - this is implicit, 1 for calc)
        # log_txt = 'Multithreading with maximum {:d} threads\n'.format(parent.threadpool.maxThreadCount())
        # parent.log.append(log_txt, alert=False)
        
        # Connect Toolbar Functions
        parent.action_Run.triggered.connect(self.start_threads)
        parent.action_Abort.triggered.connect(self.abort_workers)
        
    def start_threads(self):
        parent = self.parent
        parent.path_set.optimized_mech()
        if parent.directory.invalid: return
        if parent.optimize_running: return
        if len(parent.series_viewer.data_table) == 0: return
        if not parent.load_full_series_box.isChecked(): # TODO: May want to remove this limitation in future
            parent.log.append('"Load Full Series Into Memory" must be checked for optimization\n')
            return
        
        # Specify coefficients to be optimized
        self.coef_opt = coef_opt = self._set_coef_opt()
        if not coef_opt: return         # if nothing to optimize, don't!
                
        # Set shocks to be run
        self.shocks2run = []
        for series in parent.series.shock:
            for shock in series:
                # skip not included or exp_data not loaded from experiment
                if not shock['include'] or 'exp_data' in shock['err']: 
                    shock['SIM'] = None
                    continue
                
                # if weight variables aren't set, update
                weight_var = [shock[key] for key in ['weight_max', 'weight_min', 'weight_shift', 
                      'weight_k']]
                if np.isnan(np.hstack(weight_var)).any():
                    presize = np.shape(shock['exp_data'])[0]
                    parent.weight.update(shock=shock)
                    shock['weights'] = parent.series.weights(shock['exp_data'][:,0], shock)
                
                # if reactor temperature and pressure aren't set, update
                if np.isnan([shock['T_reactor'], shock['P_reactor']]).any():
                    parent.series.set('zone', shock['zone'])
                    
                parent.series.rate_bnds(shock)
                
                self.shocks2run.append(shock)
        
        if len(self.shocks2run) == 0: return    # if no shocks to run return
        
        shock_conditions = {'T2': [], 'P2': [], 'T5': [], 'P5': [], 'thermo_mix': []}
        for shock in self.shocks2run:
            for shock_condition in shock_conditions:
                shock_conditions[shock_condition].append(shock[shock_condition])
        
        # Set conditions of rates to be fit for each coefficient
        rxn_coef_opt = self._set_rxn_coef_opt(shock_conditions)
        
        parent.multiprocessing = parent.multiprocessing_box.isChecked()
        
        parent.update_user_settings()
        # parent.set_weights()
        
        parent.abort = False
        parent.optimize_running = True
        
        # Create mechs and duplicate mech variables
        if parent.multiprocessing == True:
            cpu_count = mp.cpu_count()
            # if cpu_count > 1: # leave open processor
                # cpu_count -= 1
            parent.max_processors = np.min([len(self.shocks2run), cpu_count])
            if parent.max_processors == 1:      # if only 1 shock, turn multiprocessing off
                parent.multiprocessing = False
            
            log_str = 'Number of processors: {:d}'.format(parent.max_processors)
            parent.log.append(log_str, alert=False)
        else:
            parent.max_processors = 1
        
        # Pass the function to execute
        self.worker = Worker(parent, self.shocks2run, parent.mech, coef_opt, rxn_coef_opt)
        self.worker.signals.result.connect(self.on_worker_done)
        self.worker.signals.finished.connect(self.thread_complete)
        self.worker.signals.update.connect(self.update)
        self.worker.signals.progress.connect(self.on_worker_progress)
        self.worker.signals.log.connect(parent.log.append)
        self.worker.signals.abort.connect(self.worker.abort)
        
        # Create Progress Bar
        # parent.create_progress_bar()
                
        if not parent.abort:
            s = 'Optimization starting\n\n   Iteration\t\t   Loss Function'
            parent.log.append(s, alert=False)
            parent.threadpool.start(self.worker)
    
    def _set_coef_opt(self):                   
        mech = self.parent.mech
        coef_opt = []
        for rxnIdx in range(mech.gas.n_reactions):      # searches all rxns
            if not mech.rate_bnds[rxnIdx]['opt']: continue        # ignore fixed reactions
            
            # check all coefficients
            for coefIdx, (coefName, coefDict) in enumerate(mech.coeffs_bnds[rxnIdx].items()):
                if coefDict['opt']:
                    coef_opt.append({'rxnIdx': rxnIdx, 'coefIdx': coefIdx, 'coefName': coefName})
        
        return coef_opt                    
    
    def _set_rxn_coef_opt(self, shock_conditions, min_T_range=1000):
        coef_opt = deepcopy(self.coef_opt)
        mech = self.parent.mech
        rxn_coef_opt = []
        for coef in coef_opt:
            if len(rxn_coef_opt) == 0 or coef['rxnIdx'] != rxn_coef_opt[-1]['rxnIdx']:
                rxn_coef_opt.append(coef)
                rxn_coef_opt[-1]['coefIdx'] = [rxn_coef_opt[-1]['coefIdx']]
                rxn_coef_opt[-1]['coefName'] = [rxn_coef_opt[-1]['coefName']]
            else:
                rxn_coef_opt[-1]['coefIdx'].append(coef['coefIdx'])
                rxn_coef_opt[-1]['coefName'].append(coef['coefName'])
        
        T_bnds = np.array([np.min(shock_conditions['T2']), np.max(shock_conditions['T2'])])
        if T_bnds[1] - T_bnds[0] < min_T_range:  # if T_range isn't large enough increase it
            T_mean = np.mean(T_bnds)
            T_bnds = np.array([T_mean-min_T_range/2, T_mean+min_T_range/2])
            # T_bnds = np.ones_like(T_bnds)*np.mean(T_bnds) + np.ones_like(T_bnds)*[-1, 1]*min_T_range/2
        invT_bnds = np.divide(10000, T_bnds)
        P_bnds = [np.min(shock_conditions['P2']), np.max(shock_conditions['P2'])]
        for rxn_coef in rxn_coef_opt:
            n_coef = len(rxn_coef['coefIdx'])
            rxn_coef['invT'] = np.linspace(*invT_bnds, n_coef)
            rxn_coef['T'] = np.divide(10000, rxn_coef['invT'])
            rxn_coef['P'] = np.linspace(*P_bnds, n_coef)
            rxn_coef['X'] = shock_conditions['thermo_mix'][0]   # TODO: IF MIXTURE COMPOSITION FOR DUMMY RATES MATTER CHANGE HERE
                      
        return rxn_coef_opt

    def update(self, result):
        loss_str = '{:.3e}'.format(result['loss']).replace('e+', 'e').replace('e-0', 'e-')
        self.parent.log.append('\t{:s} {:^5d}\t\t\t{:^s}'.format(
            result['type'][0].upper(), result['i'], loss_str), alert=False)
        self.parent.tree.update_coef_rate_from_opt(result['coef_opt'], result['x'])
        
        # if displayed shock isn't in shocks being optimized, calculate the new plot
        if result['ind_var'] is None and result['observable'] is None:
            self.parent.run_single()
        else:       # if displayed shock in list being optimized, show result
            self.parent.plot.signal.update_sim(result['ind_var'][:,0], result['observable'][:,0])
        
        self.parent.plot.opt.update(result['stat_plot'])
    
    def on_worker_progress(self, perc_completed, time_left):
        self.parent.update_progress(perc_completed, time_left)
    
    def thread_complete(self): pass
    
    def on_worker_done(self, result):
        parent = self.parent
        parent.optimize_running = False
        if result is None or len(result) == 0: return
        
        # update mech to optimized one
        if 'local' in result:
            update_mech_coef_opt(parent.mech, self.coef_opt, result['local']['x'])
        else:
            update_mech_coef_opt(parent.mech, self.coef_opt, result['global']['x'])
        
        for opt_type, res in result.items():
            total_shock_eval = (res['nfev']+1)*len(self.shocks2run)
            message = res['message'][:1].lower() + res['message'][1:]
            
            parent.log.append('\n{:s} {:s}'.format(opt_type.capitalize(), message))
            parent.log.append('\telapsed time:\t{:.2f}'.format(res['time']), alert=False)
            parent.log.append('\tloss function:\t{:.3e}'.format(res['fval']), alert=False)
            parent.log.append('\topt iters:\t\t{:.0f}'.format(res['nfev']+1), alert=False)
            parent.log.append('\tshock evals:\t{:.0f}'.format(total_shock_eval), alert=False)
            parent.log.append('\tsuccess:\t\t{:}'.format(res['success']), alert=False)
        
        parent.log.append('\n', alert=False)
        parent.save.chemkin_format(parent.mech.gas, parent.path_set.optimized_mech())
        parent.path_set.mech()  # update mech pulldown choices

    def abort_workers(self):
        if hasattr(self, 'worker'):
            self.worker.signals.abort.emit()
            self.parent.abort = True
            # self.parent.update_progress(100, '00:00:00') # This turns off the progress bar
