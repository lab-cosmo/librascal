from .coulomb_matrix import SortedCoulombMatrix

from ..lib import RepresentationManager

_representations_dict = {"coulomb":["distance","rownorm"]}


_representations = {}
for k,v in RepresentationManager.__dict__.items():
    if "pybind11_builtins.pybind11_type" in str(type(v)):
        kl = k.lower()
        for name,options in _representations_dict.items():
            for opt in options:
                if name in kl:
                    if opt in kl:
                        _representations[(name,opt)] = v

def RepresentationFactory(name,options,args):
    key = tuple([name]+sorted(options)))
    return _representations[key](*args)





