from ..lib import RepresentationManager

_representations_list = ["coulomb"]


_representations = {}
for k, v in RepresentationManager.__dict__.items():
    if "pybind11_builtins.pybind11_type" in str(type(v)):
        kl = k.lower()
        for name in _representations_list:
            if name in kl:
                _representations[name] = v


def RepresentationFactory(name, *options):
    return _representations[name]
