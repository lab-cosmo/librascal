from ..lib._rascal.models import kernels
from ..neighbourlist import AtomsList
from ..utils import return_deepcopy

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

    @return_deepcopy
    def get_init_params(self):
        init_params = dict(representation=self._rep,
            name = self.name,
            kernel_type = self.kernel_type,
            target_type = self.target_type,
            kwargs = self._kwargs)
        return init_params

    def _set_data(self, data):
        pass

    def _get_data(self):
        return dict(pseudo_points="")

    def extend(self, atoms_list, selected_indices):
        if isinstance(atoms_list, AtomsList):
            self._pseudo_points.extend(self._representation, atoms_list.managers, selected_indices)
        else:
            self._pseudo_points.extend(self._representation, atoms_list, selected_indices)
    def size(self):
        return self._pseudo_points.size()
    def get_features(self):
        return self._pseudo_points.get_features()
