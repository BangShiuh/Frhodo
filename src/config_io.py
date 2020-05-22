#!/usr/bin/env python3
# -*- coding: utf-8 -*-

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

import numpy as np

import pathlib, sys, collections

def FlowMap(*args, **kwargs):
    m = yaml.comments.CommentedMap(*args, **kwargs)
    m.fa.set_flow_style()
    return m

def FlowList(*args, **kwargs):
    lst = yaml.comments.CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return lst

# Improved float formatting requires Numpy >= 1.14
if hasattr(np, 'format_float_positional'):
    def float2string(data):
        if data == 0:
            return '0.0'
        elif 0.01 <= abs(data) < 10000:
            return np.format_float_positional(data, trim='0')
        else:
            return np.format_float_scientific(data, trim='0')
else:
    def float2string(data):
        return repr(data)

def represent_float(self, data):
    # type: (Any) -> Any
    if data != data:
        value = '.nan'
    elif data == self.inf_value:
        value = '.inf'
    elif data == -self.inf_value:
        value = '-.inf'
    else:
        value = float2string(data)

    return self.represent_scalar(u'tag:yaml.org,2002:float', value)

yaml.RoundTripRepresenter.add_representer(float, represent_float)

def deep_convert_dict(layer):   # convert all OrderedDict into dict to remove comments
    to_ret = layer              # they add a space each time prrogram is opened
    if isinstance(layer, collections.OrderedDict):
        to_ret = dict(layer)

    try:
        for key, value in to_ret.items():
            to_ret[key] = deep_convert_dict(value)
    except AttributeError:
        pass

    return to_ret

class GUI_Config(yaml.YAML):
    def __init__(self):
        super().__init__()
        self.default_flow_style = False
        self.block_seq_indent = 2
        # self.indent = 4
        self.allow_unicode = True
        self.encoding = 'utf-8'
        self.width = 80
        
        self.loader = yaml.RoundTripLoader
        
        self.setDefault()
        
    def setDefault(self):
        self.settings = {'Directory Settings': {
                            'directory file': '',
                            },
                         'Experiment Settings': {
                            'temperature units': {'zone 1': 'K',    'zone 2': 'K',    'zone 5': 'K'},
                            'pressure units':    {'zone 1': 'Torr', 'zone 2': 'Torr', 'zone 5': 'atm'},
                            'velocity units': 'm/s',
                            },
                         'Reactor Settings': {
                            'reactor': 'Incident Shock Reactor',
                            'solve energy': True,
                            'frozen composition': False,
                            'simulation end time': {'value': 12.0, 'units': 'us'},
                            'ODE solver': 'BDF',
                            'simulation interpolation factor': 1,
                            'ODE tolerance': {'relative': 1E-6, 'absolute': 1E-8},
                            },
                         'Optimization Settings': {
                            'time uncertainty': 0.0,
                            'loss function alpha': -2.00,
                            'loss function c': 1.00,
                            'multiprocessing': True,
                            'enabled':                  {'global': True,     'local': True},
                            'algorithm':                {'global': 'DIRECT', 'local': 'Subplex'},
                            'initial step':             {'global': 1.0E-2,   'local': 1.0E-2},
                            'relative tolerance x':     {'global': 1.0E-4,   'local': 1.0E-4},
                            'relative tolerance fcn':   {'global': 5.0E-4,   'local': 1.0E-3},
                            'weight function': {
                                'max': 100,
                                'min': [0, 0],
                                'time location': [0.5, 3.7],
                                'inverse growth rate': [0, 0.3],
                                },
                            },
                         'Plot Settings': {
                            'x-scale': 'linear',
                            'y-scale': 'abslog',
                            },
                        }
                         
    def to_yaml(self, dest=None):
        settings = self.settings
        out = yaml.comments.CommentedMap(settings)
        
        # reformat certain sections 
        toFlowMap = [['Experiment Settings', 'temperature units'],
                     ['Experiment Settings', 'pressure units'],
                    ]
        toFlowList = [['Optimization Settings', 'weight function', 'min'],
                      ['Optimization Settings', 'weight function', 'time location'],
                      ['Optimization Settings', 'weight function', 'inverse growth rate'],
                     ]
        for FlowType, toFlow in {'Map': toFlowMap, 'List': toFlowList}.items():
            for keys in toFlow:
                out_element = out
                settings_element = settings
                for key in keys[:-1]:
                    out_element = out_element[key]
                    settings_element = settings_element[key]
                
                if FlowType == 'Map':
                    out_element[keys[-1]] = FlowMap(settings_element[keys[-1]])
                elif FlowType == 'List':
                    out_element[keys[-1]] = FlowList(settings_element[keys[-1]])
        
        # add spacing between main sections
        for key in list(self.settings.keys())[1:]:
            out.yaml_set_comment_before_after_key(key, before='\n')
        
        # if node.note:
            # note = textwrap.dedent(node.note.rstrip())
            # if '\n' in note:
                # note = yaml.scalarstring.PreservedScalarString(note)
            # out['note'] = note

        # self.dump(representer.represent_dict(out), dest)
        if dest is None: 
            self.dump(out, sys.stdout)
        else:
            with open(dest, 'w') as configFile:
                self.dump(out, configFile)

    def from_yaml(self, src=None):
        if src is None: return
        if not src.exists(): return
        
        with open(src, 'r') as configFile:
            data = deep_convert_dict(self.load(configFile))

        self.settings.update(data)


class GUI_settings:
    def __init__(self, parent):
        self.parent = parent    # Need to find a better solution than passing parent
        self.cfg_io = GUI_Config()
        self.cfg = self.cfg_io.settings
    
    def load(self):
        parent = self.parent
        
        self.cfg_io.from_yaml(parent.path['default_config'])
        
        parent.shock_choice_box.setValue(1)
        # parent.time_offset_box.setValue(float(self.config['Experiment Settings']['time_offset']))
        # parent.time_unc_box.setValue(float(self.config['Experiment Settings']['time_unc']))
        # parent.start_ind_box.setValue(float(self.config['Experiment Settings']['start_ind']))
        # parent.weight_k_box.setValue(float(self.config['Experiment Settings']['weight_k']))
        # parent.weight_shift_box.setValue(float(self.config['Experiment Settings']['weight_shift']))
        # parent.weight_min_box.setValue(float(self.config['Experiment Settings']['weight_min'])*100)
        
        # parent.path['path_file'] = pathlib.Path(cfg['Directory Settings']['directory file'])
        parent.path_file_box.setPlainText(str(self.cfg['Directory Settings']['directory file']))
    
    def save(self, save_all=False):
        parent = self.parent
        
        self.cfg['Directory Settings']['directory file'] = str(parent.path['path_file'])
        
        '''
        if save_all:
            self.config['Directory File'] = {'file': str(parent.path['path_file'])}
            
            # self.config['Experiment Settings'] = {'start_ind':       parent.var['start ind'],
                                                  # 'time_offset':     parent.var['time_offset'],
                                                  # 'time_unc':        parent.var['time_unc'],
                                                  # 'weight_k':        parent.var['weight_k'],
                                                  # 'weight_shift':    parent.var['weight_shift'],
                                                  # 'weight_min':      parent.var['weight_min']}
                    
        else:
            self.config.set('Directory File', 'file', str(parent.path['path_file']))
        '''
        
        self.cfg_io.to_yaml(parent.path['default_config'])
        

if __name__ == '__main__':
    gui_cfg = GUI_Config()

    path = {}
    path['default_config'] = pathlib.Path('default_config.yaml')

    # gui_cfg.dump(gui_cfg.settings, sys.stdout)
    gui_cfg.to_yaml(path['default_config'])
    gui_cfg.from_yaml(path['default_config'])