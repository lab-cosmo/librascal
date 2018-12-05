from ..lib import RepresentationManager

_representations_dict = {
    "coulomb":[["distance","rownorm"]],
    "sphericalexpansion":[["expansion"]]
    }


_representations = {}
for k,v in RepresentationManager.__dict__.items():
    if "pybind11_builtins.pybind11_type" in str(type(v)):
        kl = k.lower()
        for name,options_list in _representations_dict.items():
            for options in options_list:
                for opt in options:
                  if name in kl:
                      if opt in kl:
                          _representations[(name,opt)] = v

def RepresentationFactory(name,*options):
    key = tuple([name]+list(options))
    return _representations[key]