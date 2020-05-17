#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file is part of Frhodo. Copyright © 2020, UChicago Argonne, LLC
# and licensed under BSD-3-Clause. See License.txt in the top-level 
# directory for license and copyright information.

import os, sys, platform, multiprocessing, pathlib, logging, traceback
from logging.handlers import RotatingFileHandler
# os.environ['QT_API'] = 'pyside2'        # forces pyside2

from qtpy.QtWidgets import QMainWindow, QApplication, QWidget, QDialog
from qtpy import uic, QtCore, QtGui

import numpy as np
# from timeit import default_timer as timer

import plot, misc_widget, options_panel_widgets, convert_units, sim_explorer_widget
import mech_fcns, settings, save_widget
    
if os.environ['QT_API'] == 'pyside2': # Silence warning: "Qt WebEngine seems to be initialized from a plugin."
    QApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    
# Handle high resolution displays:  Minimum recommended resolution 1280 x 960
if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

class Main(QMainWindow):
    def __init__(self):
        super().__init__()
        self.path_set = settings.path(self)
        uic.loadUi(str(self.path['main']/'UI'/'main_window.ui'), self)
        self.splitter.moveSplitter(0, 1)    # moves splitter 0 as close to 1 as possible
        self.setWindowIcon(QtGui.QIcon(str(self.path['main']/'UI'/'graphics'/'main_icon.png')))
        
        # Start threadpools
        self.threadpool = QtCore.QThreadPool()
        self.threadpool.setMaxThreadCount(2) # Sets thread count to 1 (1 for gui - this is implicit, 1 for calc)
        
        # Set selected tabs
        for tab_widget in [self.option_tab_widget, self.plot_tab_widget]:
            tab_widget.setCurrentIndex(0)
        
        # Set Clipboard
        self.clipboard = QApplication.clipboard()
        
        self.var = {'reactor': {'t_unit_conv': 1}}
        self.SIM = mech_fcns.Simulation_Result()
        self.mech_loaded = False
        self.convert_units = convert_units.Convert_Units(self)
        self.series = settings.series(self)
        
        self.sim_explorer = sim_explorer_widget.SIM_Explorer_Widgets(self)
        self.plot = plot.All_Plots(self)
        options_panel_widgets.Initialize(self)
        self.mech = mech_fcns.Chemical_Mechanism()
        
        # Setup save sim
        self.save_sim = save_widget.Save_Dialog(self)
        self.save_sim_button.clicked.connect(self.save_sim.execute)
        self.action_Save.triggered.connect(self.save_sim.execute)
        
        # Initialize Settings
        self.initialize_settings()
        
        self.show()
    
    def initialize_settings(self):
        self.var['old_shock_choice'] = self.var['shock_choice'] = 1
                
        self.user_settings = settings.user_settings(self)
        if self.path['default_settings.ini'].exists():
            self.user_settings.load(self.path['default_settings.ini'])
        
        self.load_full_series = self.load_full_series_box.isChecked()   # TODO: Move to somewhere else?
        
        # load previous paths if file in path, can be accessed, and is a file
        if ('path_file' in self.path and os.access(self.path['path_file'], os.R_OK) and 
            self.path['path_file'].is_file()):
            
            self.path_set.load_dir_file(self.path['path_file'])
            
        self.update_user_settings()
    
    def load_mech(self, event = None):
        def mechhasthermo(mech_path):
            f = open(mech_path, 'r')
            while True:
                line = f.readline()
                if '!' in line[0:2]:
                    continue
                if 'thermo' in line.split('!')[0].strip().lower():
                    return True
                
                if not line:
                    break
                
            f.close()
            return False
        
        if self.mech_select_comboBox.count() == 0: return   # if no items return, unsure if this is needed now
        
        # Specify mech file path
        self.path['mech'] = self.path['mech_main'] / str(self.mech_select_comboBox.currentText())
        if not self.path['mech'].is_file(): # if it's not a file, then it was deleted
            self.path_set.mech()            # update mech pulldown choices
            return
        
        # Check use thermo box viability
        if mechhasthermo(self.path['mech']):
            self.use_thermo_file_box.setEnabled(True)
            # Autoselect checkbox off if thermo exists in mech
            if self.sender() is None or 'use_thermo_file_box' not in self.sender().objectName(): 
                self.use_thermo_file_box.setChecked(False)      
        else:
            self.use_thermo_file_box.setChecked(True)
            self.use_thermo_file_box.setDisabled(True) # disable checkbox if no thermo in mech file
        
        # Enable thermo select based on use_thermo_file_box
        if self.use_thermo_file_box.isChecked():
            self.thermo_select_comboBox.setEnabled(True)
        else:
            self.thermo_select_comboBox.setDisabled(True)
        
        # Specify thermo file path        
        if self.use_thermo_file_box.isChecked():
            if self.thermo_select_comboBox.count() > 0:
                self.path['thermo'] = self.path['mech_main'] / str(self.thermo_select_comboBox.currentText())
            else:
                self.log.append('Error loading mech:\nNo thermodynamics given', alert=True)
                return
        else:
            self.path['thermo'] = None
        
        # Initialize Mechanism
        self.log.clear([])  # Clear log when mechanism changes to avoid log errors about prior mech
        mech_load_output = self.mech.load_mechanism(self.path)
        self.log.append(mech_load_output['message'], alert=not mech_load_output['success'])
        self.mech_loaded = mech_load_output['success']
        
        if not mech_load_output['success']:   # if error: update species and return
            self.mix.update_species()
            self.log._blink(True)   # updating_species is causing blink to stop due to successful shock calculation
            return
        
        # Initialize tables and trees
        self.tree.set_trees(self.mech)
        self.mix.update_species()       # this was commented out, could be because multiple calls to solver from update_mix / setItems
        
        tabIdx = self.plot_tab_widget.currentIndex()
        tabText = self.plot_tab_widget.tabText(tabIdx)
        if tabText == 'Signal/Sim':
            # Force observable_widget to update
            observable = self.plot.observable_widget.widget['main_parameter'].currentText()
            self.plot.observable_widget.widget['main_parameter'].currentIndexChanged[str].emit(observable)
        elif tabText == 'Sim Explorer': # TODO: This gets called twice?
            self.sim_explorer.update_all_main_parameter()
            
    def shock_choice_changed(self, event):
        if 'exp_main' in self.directory.invalid:    # don't allow shock change if problem with exp directory
            return
            
        self.var['old_shock_choice'] = self.var['shock_choice']
        self.var['shock_choice'] = event
        
        self.shockRollingList = ['P1', 'u1']    # reset rolling list
        self.rxn_change_history = []  # reset tracking of rxn numbers changed
            
        if not self.optimize_running:
            self.log.clear([])
        self.series.change_shock()  # link display_shock to correct set and 
    
    def update_user_settings(self, event = None):
        # This is one is located on the Files tab
        shock = self.display_shock
        self.series.set('series_name', self.exp_series_name_box.text())
        
        t_unit_conv = self.var['reactor']['t_unit_conv']
        if self.time_offset_box.value()*t_unit_conv != shock['time_offset']: # if values are different
            self.series.set('time_offset', self.time_offset_box.value()*t_unit_conv)
            if hasattr(self.mech_tree, 'rxn'):              # checked if copy valid in function
                self.tree._copy_expanded_tab_rates()        # copy rates and time offset
        
        self.var['time_unc'] = self.time_unc_box.value()*t_unit_conv
                
        if event is not None:
            sender = self.sender().objectName()
            if 'time_offset' in sender and hasattr(self, 'SIM'): # Don't rerun SIM if it exists
                if hasattr(self.SIM, 'independent_var') and hasattr(self.SIM, 'observable'):
                    self.plot.signal.update_sim(self.SIM.independent_var, self.SIM.observable)
            elif any(x in sender for x in ['end_time', 'sim_interp_factor', 'ODE_solver', 'rtol', 'atol']):
                self.run_single()
            elif self.display_shock['exp_data'].size > 0: # If exp_data exists
                self.plot.signal.update(update_lim=False)
                self.plot.signal.canvas.draw()
        '''
        # debug
        for i in self.var:
            print('key: {:<14s} value: {:<16s}'.format(i, str(self.var[i])))
        '''
    
    def keyPressEvent(self, event): pass
        # THIS IS NOT FULLY FUNCTIONING
        # http://ftp.ics.uci.edu/pub/centos0/ics-custom-build/BUILD/PyQt-x11-gpl-4.7.2/doc/html/qkeyevent.html
        # print(event.modifiers(),event.text())
    
    def run_single(self, event=None, t_save=None, rxn_changed=False):
        if not self.mech_loaded: return                 # if mech isn't loaded successfully, exit
        if not hasattr(self.mech_tree, 'rxn'): return   # if mech tree not set up, exit
        
        shock = self.display_shock
        
        T_reac, P_reac, mix = shock['T_reactor'], shock['P_reactor'], shock['thermo_mix']
        self.tree.update_rates()
        
        # calculate all properties or observable by sending t_save
        tabIdx = self.plot_tab_widget.currentIndex()
        tabText = self.plot_tab_widget.tabText(tabIdx)
        if tabText == 'Sim Explorer':
            t_save = np.array([0])
        
        SIM_kwargs = {'u_reac': shock['u2'], 'rho1': shock['rho1'], 'observable': self.display_shock['observable'], 
            't_lab_save': t_save, 'sim_int_f': self.var['reactor']['sim_interp_factor'], 
            'ODE_solver': self.var['reactor']['ode_solver'], 
            'rtol': self.var['reactor']['ode_rtol'], 'atol': self.var['reactor']['ode_atol']}
        
        if '0d Reactor' in self.var['reactor']['name']:
            SIM_kwargs['solve_energy'] = self.var['reactor']['solve_energy']
            SIM_kwargs['frozen_comp'] = self.var['reactor']['frozen_comp']
        
        self.SIM, verbose = self.mech.run(self.var['reactor']['name'], self.var['reactor']['t_end'], 
                                          T_reac, P_reac, mix, **SIM_kwargs)
        
        if verbose['success']:
            self.log._blink(False)
        else:
            self.log.append(verbose['message'])
        
        if self.SIM is not None:
            self.plot.signal.update_sim(self.SIM.independent_var, self.SIM.observable, rxn_changed)
            if tabText == 'Sim Explorer':
                self.sim_explorer.populate_main_parameters()
                self.sim_explorer.update_plot(self.SIM) # somtimes duplicate updates
        else:
            nan = np.array([np.nan, np.nan])
            self.plot.signal.update_sim(nan, nan)   # make sim plot blank
            if tabText == 'Sim Explorer':
                self.sim_explorer.update_plot(None)
            return # If mech error exit function

    # def raise_error(self):
        # assert False

class Error_Window(QDialog):
    def __init__(self, path):
        super().__init__()
        self.path = path
        uic.loadUi(str(self.path['main']/'UI'/'error_window.ui'), self)
        self.setWindowIcon(QtGui.QIcon(str(self.path['main']/'UI'/'graphics'/'main_icon.png')))
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowTitleHint |
                            QtCore.Qt.WindowCloseButtonHint | QtCore.Qt.WindowStaysOnTopHint)

        self.close_button.clicked.connect(self.closeEvent) 
        self.installEventFilter(self)
        
        self.exec_()
        
    def eventFilter(self, obj, event):
        # intercept enter, space and escape
        if event.type() == QtCore.QEvent.KeyPress:
            if event.key() in [QtCore.Qt.Key_Escape, QtCore.Qt.Key_Return, QtCore.Qt.Key_Space]:
                self.close_button.click()
                return True
                    
        return super().eventFilter(obj, event)
    
    def closeEvent(self, event):
        QApplication.quit() # some errors can be recovered from, maybe I shouldn't autoclose the program

def excepthook(type, value, tback):
    # log the exception
    path = {'main': pathlib.Path(sys.argv[0]).parents[0].resolve()}
    path['log'] = path['main']/'error.log'
    
    log_formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    
    log_handler = RotatingFileHandler(path['log'], mode='a', maxBytes=1*1024*1024, # maximum of 1 MB
                                     backupCount=1, encoding=None, delay=0)        # maximum of 2 error files
    log_handler.setFormatter(log_formatter)
    log_handler.setLevel(logging.DEBUG)

    app_log = logging.getLogger('root')
    app_log.setLevel(logging.DEBUG)
    app_log.addHandler(log_handler)
    
    text = "".join(traceback.format_exception(type, value, tback))   
    app_log.error(text)

    # call the default handler
    sys.__excepthook__(type, value, tback)
    
    Error_Window(path)    

sys.excepthook = excepthook
            
if __name__ == '__main__':
    if platform.system() == 'Windows':  # this is required for pyinstaller on windows
        multiprocessing.freeze_support()

    app = QApplication(sys.argv)
    main = Main()
    sys.exit(app.exec_())
   