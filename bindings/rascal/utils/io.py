import os
import importlib
from collections import Iterable
import numpy as np
import json

BETA_VERSION = '0.1'

CURRENT_VERSION = BETA_VERSION

def dump_obj(fn,instance,version=CURRENT_VERSION):
    if isinstance(instance, BaseIO):
        instance.to_file(fn,version)
    else:
        raise RuntimeError('The instance does not inherit from BaseIO: {}'.format(instance.__class__.__mro__))

def load_obj(fn):
    return BaseIO().from_file(fn)

def from_dict(data):
    return BaseIO().from_dict(data)

def dump_json(fn,data):
    with open(fn,'w') as f:
        json.dump(data,f,sort_keys=True,indent=2)

def load_json(fn):
    with open(fn,'r') as f:
        data = json.load(f)
    return data

def is_npy(data):
    return isinstance(data, np.ndarray)

def is_large_array(data):
    if is_npy(data):
        if data.nbytes > 50e6:
            return True
        else:
            return False
    else:
        return False

def is_npy_filename(fn):
    if isinstance(fn,str):
        filename , file_extension = os.path.splitext(fn)
        if file_extension == '.npy':
            return True
        else:
            return False
    else:
        return False

def get_class(module_name, class_name):
    module = importlib.import_module(module_name)
    class_ = getattr(module, class_name)
    return class_

def obj2dict_beta(cls,state):
    VERSION = BETA_VERSION
    module_name = cls.__module__
    class_name = cls.__name__
    frozen = dict(version=VERSION, class_name=class_name,
                    module_name=module_name,
                    init_params=state['init_params'],
                    data=state['data'])
    return frozen

def dict2obj_beta(data):
    cls = get_class(data['module_name'], data['class_name'])
    obj = cls(**data['init_params'])
    obj._set_data(data['data'])
    return obj

def is_valid_object_dict_beta(data):
    valid_keys = ['version','class_name','module_name','init_params','data']
    aa = []
    if isinstance(data,dict):
        for k in data:
            if k in valid_keys:
                aa.append(True)
        if len(aa) == len(valid_keys):
            return True
        else:
            return False
    else:
        return False

# def is_a_wrapped_object(obj):
#     """Answer the question: is obj imported from C++ binding library ?
#     namely from 'rascal.lib._rascal'."""
#     module_str = getattr(obj, '__module__', None)
#     module_root_str = '.'.join(module_str.split('.')[:3])
#     return 'rascal.lib._rascal' == module_root_str

obj2dict = {BETA_VERSION:obj2dict_beta}
dict2obj = {BETA_VERSION:dict2obj_beta}
is_valid_object_dict = {BETA_VERSION:is_valid_object_dict_beta}

class BaseIO(object):
    """
    Base class giving arbitrary python class the ability to serialize themselves
    into python dictionary and being saved to files.
    A class that inherits from BaseIO should implement the _get_data and get_init_params functions (see _get_state() for more details).
    """
    def __init__(self):
        super(BaseIO,self).__init__()

    def _get_state(self):
        """fetch the state of self using the _get_data and get_init_params functions that should be implemented in the Implementation class that inherits from BaseIO.

        get_init_params is expected to return a dictionary containing all the
        parameters used by the __init__() methods.

        _get_data is expected to return a dictionary containing all the data
        that is not set by the initialization of the class."""
        state = dict(data=self._get_data(),
                     init_params=self.get_init_params())
        return state

    def to_dict(self,obj=None,version=CURRENT_VERSION):
        """recursirvely convert the python object obj into a dictionary that
        can be used to create a copy of obj.
        obj has to inherit from BaseIO and implements a _get_data and
        get_init_params methods"""
        if obj is None:
            obj = self
        state = obj._get_state()

        for name,entry in state.items():
            if isinstance(entry, dict):
                for k,v in entry.items():
                    if isinstance(v, BaseIO) is True:
                        state[name][k] = self.to_dict(v,version)
                    elif isinstance(v, list):
                        ll = []
                        for val in v:
                            if isinstance(val, BaseIO) is True:
                                ll.append(self.to_dict(val,version))
                            else:
                                ll.append(val)
                        state[name][k] = ll

        data = obj2dict[version](obj.__class__, state)
        return data

    def from_dict(self,data):
        """recursirvely convert python a dictionary that
        obj has to inherit from BaseIO."""
        version = data['version']
        for name,entry in data.items():
            if isinstance(entry, dict):
                for k,v in entry.items():
                    if is_valid_object_dict[version](v) is True:
                         data[name][k] = self.from_dict(v)
                    elif isinstance(v, list):
                        ll = []
                        for val in v:
                            if is_valid_object_dict[version](val) is True:
                                ll.append(self.from_dict(val))
                            else:
                                ll.append(val)

                        data[name][k] = ll


        obj = dict2obj[version](data)
        return obj

    def to_file(self,fn,version=CURRENT_VERSION):
        fn = os.path.abspath(fn)
        filename , file_extension = os.path.splitext(fn)
        data = self.to_dict(version=version)
        class_name = data['class_name'].lower()

        if file_extension == '.json':
            self._dump_npy(fn,data,class_name)
            dump_json(fn,data)
        else:
            raise NotImplementedError('Unknown file extention: {}'.format(file_extension))

    def from_file(self,fn):
        fn = os.path.abspath(fn)
        path = os.path.dirname(fn)
        filename , file_extension = os.path.splitext(fn)
        if file_extension == '.json':
            data = load_json(fn)
            version = data['version']
            if is_valid_object_dict[version](data):
                self._load_npy(data, path)
                return self.from_dict(data)
            else:
                raise RuntimeError('The file: {}; does not contain a valid dictionary representation of an object.'.format(fn))

        else:
            raise NotImplementedError('Unknown file extention: {}'.format(file_extension))

    def _dump_npy(self,fn,data,class_name):
        filename , file_extension = os.path.splitext(fn)
        for k,v in data.items():
            if isinstance(v,dict):
                if 'class_name' in data:
                    class_name = data['class_name'].lower()
                self._dump_npy(fn, v, class_name)
            elif is_large_array(v) is True:
                if 'tag' in data:
                    class_name += '-'+ data['tag']
                v_fn = filename + '-{}-{}'.format(class_name,k) + '.npy'
                v_bfn = os.path.basename(v_fn)
                data[k] = v_bfn
                np.save(v_fn,v)

            elif is_npy(v) is True:
                data[k] = ['npy',v.tolist()]

    def _load_npy(self,data, path):
        for k,v in data.items():
            if isinstance(v,dict):
                self._load_npy(v, path)
            elif is_npy_filename(v) is True:
                data[k] = np.load(os.path.join(path, v),mmap_mode='r')
            elif isinstance(v,list):
                if len(v) == 2:
                    if 'npy' == v[0]:
                        data[k] = np.array(v[1])

