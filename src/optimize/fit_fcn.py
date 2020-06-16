# This file is part of Frhodo. Copyright © 2020, UChicago Argonne, LLC
# and licensed under BSD-3-Clause. See License.txt in the top-level 
# directory for license and copyright information.

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.interpolate import interp1d
from scipy import stats
from copy import deepcopy

import mech_fcns
from optimize.fit_coeffs import fit_coeffs

mpMech = {}

def initialize_parallel_worker(mech_txt, coeffs, coeffs_bnds, rate_bnds):
    mpMech['obj'] = mech = mech_fcns.Chemical_Mechanism()
    mech.set_mechanism(mech_txt)    # load mechanism from yaml text in memory
    mech.coeffs = deepcopy(coeffs)
    mech.coeffs_bnds = deepcopy(coeffs_bnds)
    mech.rate_bnds = deepcopy(rate_bnds)

def outlier(res, a=2, c=1, weights=[], max_iter=25):
    def diff(res_outlier):
        if len(res_outlier) < 2: 
            return 1
        else:
            return np.diff(res_outlier)[0]
    
    if a != 2: # define outlier with 1.5 IQR rule
        trunc_res = np.abs(res.copy())
        percentile = 25
        if len(weights) == len(res):
            trunc_res = trunc_res[weights > 0.95*np.max(weights)]  # This computes the outlier threshold based on weights >= 0.95
                
        res_outlier = []
        for n in range(max_iter):
            if diff(res_outlier) == 0:   # iterate until res_outlier is the same as prior iteration
                break
                
            if len(res_outlier) > 0:
                trunc_res = trunc_res[trunc_res < res_outlier[-1]] 
                
            q1, q3 = np.nanpercentile(trunc_res, percentile), np.nanpercentile(trunc_res, 100-percentile)
            iqr = q3 - q1       # interquartile range      
            
            if len(res_outlier) == 2:
                del res_outlier[0]
            
            res_outlier.append(q3 + iqr*1.5)
        
        res_outlier = res_outlier[-1]
                
    else:
        res_outlier = 1
        
    return c*res_outlier
    
def generalized_loss_fcn(res, a=2, c=1, weights=[]):    # defaults to L2 loss
    x_c_2 = np.power(res/c, 2)
    if a == 2:
        loss = 0.5*x_c_2
    elif a == 0:
        loss = np.log(0.5*x_c_2+1)
    elif a <= -1000:  # supposed to be negative infinity
        loss = 1 - np.exp(-0.5*x_c_2)
    else:
        loss = np.abs(a-2)/a*(np.power(x_c_2/np.abs(a-2) + 1, a/2) - 1)
        
    if len(weights) == len(loss):
        loss = np.multiply(loss, weights)
        
    return loss*np.abs(c)   # multiplying by c is not necessary, but makes order appropriate

def update_mech_coef_opt(mech, coef_opt, x):
    mech_changed = False
    for i, idxDict in enumerate(coef_opt):
        rxnIdx, coefName = idxDict['rxnIdx'], idxDict['coefName']
        if mech.coeffs[rxnIdx][coefName] != x[i]:       # limits mech changes. Should increase speed a little
            mech_changed = True
            mech.coeffs[rxnIdx][coefName] = x[i]
    
    if mech_changed:
        mech.modify_reactions(mech.coeffs)  # Update mechanism with new coefficients
  
def calculate_residuals(args_list):   
    def calc_exp_bounds(t_sim, t_exp):
        t_bounds = [np.max([t_sim[0], t_exp[0]])]       # Largest initial time in SIM and Exp
        t_bounds.append(np.min([t_sim[-1], t_exp[-1]])) # Smallest final time in SIM and Exp
        # Values within t_bounds
        exp_bounds = np.where(np.logical_and((t_exp >= t_bounds[0]),(t_exp <= t_bounds[1])))[0]
        
        return exp_bounds
    
    def time_adjust_func(t_offset, t_adjust, t_sim, obs_sim, t_exp, obs_exp, weights, verbose=False):
        t_sim_shifted = t_sim + t_offset + t_adjust

        # Compare SIM Density Grad vs. Experimental
        exp_bounds = calc_exp_bounds(t_sim_shifted, t_exp)
        t_exp, obs_exp, weights = t_exp[exp_bounds], obs_exp[exp_bounds], weights[exp_bounds]
        
        f_interp = interp1d(t_sim_shifted.flatten(), obs_sim.flatten(), kind = 'cubic')
        obs_sim_interp = f_interp(t_exp)
        
        resid = np.subtract(obs_exp, obs_sim_interp)
        resid_outlier = outlier(resid, a=var['loss_alpha'], c=var['loss_c'], 
                                weights=weights)
        
        if verbose:
            output = {'resid': resid, 'resid_outlier': resid_outlier,
                      'weights': weights,
                      'obs_sim_interp': obs_sim_interp}
            
            return output
        else:   # needs to return single value for optimization
            return generalized_loss_fcn(resid, a=var['loss_alpha'], c=resid_outlier, 
                                        weights=weights).sum()
    
    def calc_density(x, data, dim=1):
        stdev = np.std(data)
        A = np.min([np.std(data), stats.iqr(data)/1.34])/stdev  # bandwidth is multiplied by std of sample
        bw = 0.9*A*len(data)**(-1./(dim+4))

        return stats.gaussian_kde(data, bw_method=bw)(x)
        
    def OoM(x):
        if not isinstance(x, np.ndarray):
            x = np.array(x)
        x[x==0] = 1                       # if zero, make OoM 0
        return np.floor(np.log10(np.abs(x)))
    
    var, coef_opt, x, shock = args_list
    mech = mpMech['obj']
    
    # Optimization Begins, update mechanism
    update_mech_coef_opt(mech, coef_opt, x)
    
    T_reac, P_reac, mix = shock['T_reactor'], shock['P_reactor'], shock['thermo_mix']
    
    SIM_kwargs = {'u_reac': shock['u2'], 'rho1': shock['rho1'], 'observable': shock['observable'], 
                  't_lab_save': None, 'sim_int_f': var['sim_interp_factor'], 
                  'ODE_solver': var['ode_solver'], 'rtol': var['ode_rtol'], 'atol': var['ode_atol']}
    
    if '0d Reactor' in var['name']:
        SIM_kwargs['solve_energy'] = var['solve_energy']
        SIM_kwargs['frozen_comp'] = var['frozen_comp']
    
    SIM, verbose = mech.run(var['name'], var['t_end'], T_reac, P_reac, mix, **SIM_kwargs)
        
    if SIM.success:
        shock['SIM'] = SIM
    else:
        shock['SIM'] = None
    
    ind_var, obs = SIM.independent_var[:,None], SIM.observable[:,None]
    
    weights = shock['norm_weights_trim']
    obs_exp = shock['exp_data_trim']
    
    if not np.any(var['t_unc']):
        t_unc = 0
    else:
        t_unc_OoM = np.mean(OoM(var['t_unc']))  # Do at higher level? (computationally efficient)
        time_adj_decorator = lambda t_adjust: time_adjust_func(shock['time_offset'], t_adjust*10**t_unc_OoM, 
                ind_var, obs, obs_exp[:,0], obs_exp[:,1], weights)
        
        res = minimize_scalar(time_adj_decorator, bounds=var['t_unc']/10**t_unc_OoM, method='bounded')
        t_unc = res.x*10**t_unc_OoM
    
    output = time_adjust_func(shock['time_offset'], t_unc, ind_var, obs, 
                obs_exp[:,0], obs_exp[:,1], weights, verbose=True)  
    
    output['shock'] = shock
    
    plot_stats = True
    if plot_stats:
        x = np.linspace(output['resid'].min(), output['resid'].max(), 300)
        density = calc_density(x, output['resid'], dim=1)   #kernel density estimation
        output['KDE'] = np.column_stack((x, density))

    return output

# Using optimization vs least squares curve fit because y_range's change if time_offset != 0
class Fit_Fun:
    def __init__(self, input_dict):
        self.parent = input_dict['parent']
        self.shocks2run = input_dict['shocks2run']
        self.data = self.parent.series.shock
        self.coef_opt = input_dict['coef_opt']
        self.rxn_coef_opt = input_dict['rxn_coef_opt']
        self.x0 = input_dict['x0']
        self.mech = input_dict['mech']
        self.var = self.parent.var
        self.t_unc = (-self.var['time_unc'], self.var['time_unc'])
        
        self.opt_type = 'local' # this is updated outside of the class
        
        self.loss_alpha = self.parent.optimization_settings.get('loss', 'alpha')
        self.loss_c = self.parent.optimization_settings.get('loss', 'c')
        
        if 'multiprocessing' in input_dict:
            self.multiprocessing = input_dict['multiprocessing']
        
        if 'pool' in input_dict:
            self.pool = input_dict['pool']
        else:
            self.multiprocessing = False
        
        self.signals = input_dict['signals']
        
        self.i = 0        
        self.__abort = False
    
    def __call__(self, s, optimizing=True):
        def append_output(output_dict, calc_resid_output):
            for key in calc_resid_output:
                if key not in output_dict:
                    output_dict[key] = []
                    
                output_dict[key].append(calc_resid_output[key])
            
            return output_dict
        
        if self.__abort: 
            raise Exception('Optimization terminated by user')
            self.signals.log.emit('\nOptimization aborted')
            return
        
        # Convert to mech values
        x = self.fit_all_coeffs(np.exp(s*self.x0))
        if x is None: 
            return np.inf
        
        # Run Simulations
        output_dict = {}
        
        var_dict = {key: val for key, val in self.var['reactor'].items()}
        var_dict['t_unc'] = self.t_unc
        var_dict.update({'loss_alpha': 2, 'loss_c': 1}) # loss function here is for finding t_unc, mse seems to work best.
        # var_dict.update({'loss_alpha': self.loss_alpha, 'loss_c': self.loss_c})
        
        if self.multiprocessing:
            args_list = ((var_dict, self.coef_opt, x, shock) for shock in self.shocks2run)
            calc_resid_outputs = self.pool.map(calculate_residuals, args_list)
            for calc_resid_output, shock in zip(calc_resid_outputs, self.shocks2run):
                shock['SIM'] = calc_resid_output['shock']['SIM']
                append_output(output_dict, calc_resid_output)

        else:
            mpMech['obj'] = self.mech
            
            for shock in self.shocks2run:
                args_list = (var_dict, self.coef_opt, x, shock)
                calc_resid_output = calculate_residuals(args_list)
                shock['SIM'] = calc_resid_output['shock']['SIM']
                append_output(output_dict, calc_resid_output)
        
        allResid = np.concatenate(output_dict['resid'], axis=0)
        weights = np.concatenate(output_dict['weights'], axis=0)
        resid_outlier = outlier(allResid, a=self.loss_alpha, c=self.loss_c, 
                                weights=weights)
        total_loss = generalized_loss_fcn(allResid, a=self.loss_alpha, c=resid_outlier, 
                                          weights=weights).sum()    
        # For updating
        self.i += 1
        if not optimizing or self.i % 1 == 0:#5 == 0: # updates plot every 5
            if total_loss == 0:
                total_loss = np.inf
                
            shock = self.parent.display_shock
            if shock['include']:
                ind_var, observable = shock['SIM'].independent_var[:,None], shock['SIM'].observable[:,None]
            else:
                ind_var, observable = None, None
            
            stat_plot = {'shocks2run': self.shocks2run, 'resid': output_dict['resid'], 
                        'resid_outlier': resid_outlier, 'weights': output_dict['weights']}
            
            if 'KDE' in output_dict:
                stat_plot['KDE'] = output_dict['KDE']
                
                stat_plot['fit_result'] = fitres = stats.gennorm.fit(allResid)
                stat_plot['QQ'] = []
                for resid in stat_plot['resid']:
                    QQ = stats.probplot(resid, sparams=fitres, dist=stats.gennorm, fit=False)
                    QQ = np.array(QQ).T
                    stat_plot['QQ'].append(QQ)
            
            update = {'type': self.opt_type, 'i': self.i, 
                      'loss': total_loss, 'stat_plot': stat_plot, 
                      'x': x, 'coef_opt': self.coef_opt, 'ind_var': ind_var, 'observable': observable}
            
            self.signals.update.emit(update)
                  
        if optimizing:
            return total_loss
        else:
            return total_loss, x, output_dict['shock']
            
    def fit_all_coeffs(self, all_rates):      
        coeffs = []
        i = 0
        for rxn_coef in self.rxn_coef_opt:
            rxnIdx = rxn_coef['rxnIdx']
            T, P, X = rxn_coef['T'], rxn_coef['P'], rxn_coef['X']
            rxn_rates = all_rates[i:i+len(T)]
            if len(coeffs) == 0:
                coeffs = fit_coeffs(rxn_rates, T, P, X, rxn_coef['coefName'], rxnIdx, self.mech)
                if coeffs is None:
                    return
            else:
                coeffs_append = fit_coeffs(rxn_rates, T, P, X, rxn_coef['coefName'], rxnIdx, self.mech)
                if coeffs_append is None:
                    return
                coeffs = np.append(coeffs, coeffs_append)
            
            i += len(T)
        
        return coeffs