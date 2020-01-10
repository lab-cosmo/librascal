from ..lib._rascal.models import kernels
from ..neighbourlist import AtomsList

# names of existing pseudo points implementation on the pybinding side.
_pseudo_points = {}
for k, v in kernels.__dict__.items():
    if "PseudoPoints" in k:
        name = k
        _pseudo_points[name] = v

class PseudoPoints(object):
    def __init__(self, representation):
        self._representation = representation._representation
        if 'SphericalInvariants' in str(representation):
            self._pseudo_points = _pseudo_points['PseudoPointsBlockSparse_SphericalInvariants']()
        else:
            raise ValueError('No pseudo point is appropiate for '+str(representation))
    def extend(self, atoms_list, selected_indicies):
        if isinstance(atoms_list, AtomsList):
            self._pseudo_points.extend(self._representation, atoms_list.managers, selected_indicies)
        else:
            self._pseudo_points.extend(self._representation, atoms_list, selected_indicies)
    def size(self):
        return self._pseudo_points.size()
    def get_features(self):
        return self._pseudo_points.get_features()
