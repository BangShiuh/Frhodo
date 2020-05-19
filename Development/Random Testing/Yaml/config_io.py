#!/usr/bin/env python3
# -*- coding: utf-8 -*-

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

import numpy as np

import pathlib, sys

BlockMap = yaml.comments.CommentedMap

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

class GUI_Config(yaml.YAML):
    def __init__(self):
        super().__init__()
        self.default_flow_style = False
        self.block_seq_indent = 2
        # self.indent = 4
        self.allow_unicode = True
        self.encoding = 'utf-8'
        
        self.data = self.load('')
      
    def to_yaml(self, data, dest):
        out = BlockMap([('model', 'NASA7')])
        out['temperature-ranges'] = FlowList([node.Tmin, node.Tmid, node.Tmax])
        out['data'] = [FlowList(node.low_coeffs), FlowList(node.high_coeffs)]
        if node.note:
            note = textwrap.dedent(node.note.rstrip())
            if '\n' in note:
                note = yaml.scalarstring.PreservedScalarString(note)
            out['note'] = note

        self.dump(representer.represent_dict(out), dest)
        
    def prepare(self, input):
        def is_container(obj):
            if isinstance(obj, str):
                return False
            else:
                return hasattr(type(obj), '__iter__')
            
        if isinstance(input, dict):
            for k, v in input.items():
                if is_container(v):
                    print("{0} : ".format(k))
                    self.prepare(v)
                else:
                    print("{0} : {1}".format(k, v))
                    
        elif is_container(input) and not isinstance(input, dict):
            has_any_subcontainers = any(is_container(obj) for obj in input)
            if not has_any_subcontainers:
                print("{0}".format(input))
            else:
                for v in input:
                    if is_container(v):
                        self.prepare(v)
                    else:
                        print("{0}".format(v))
            

gui_cfg = GUI_Config()

path = {}
path['default_config'] = 'default_config.yaml'

settings = {'top level': {'Path': 'muckity', 'opt': {'min': [{'hi': 'nice'}, 1], 'max': [1, [11,13]]}}}

gui_cfg.prepare(settings)

# gui_cfg.dump(settings, sys.stdout)

# with open(path['default_config'], 'w') as configFile:
    # data = gui_cfg.dump(settings, configFile)
    