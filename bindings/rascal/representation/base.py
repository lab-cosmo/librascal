from ..lib import RepresentationManager

_representations_list = ["sortedcoulomb", "sphericalexpansion"]


_representations = {}
for k, v in RepresentationManager.__dict__.items():
    if "pybind11_builtins.pybind11_type" in str(type(v)):
        kl = k.lower()
        for name in _representations_list:
            if name in kl:
                _representations[kl] = v


def RepresentationFactory(name, *args):
    if name not in _representations:
        raise NameError('The representations factory {} has not been registered. The available combinations are: {}'.format(
            name, list(_representations.keys())))
    return _representations[name](*args)
