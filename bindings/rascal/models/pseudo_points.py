from ..lib._rascal.models import kernels
from ..neighbourlist import AtomsList
from ..utils import BaseIO

# names of existing pseudo points implementation on the pybinding side.
_pseudo_points = {}
for k, v in kernels.__dict__.items():
    if "PseudoPoints" in k:
        name = k
        _pseudo_points[name] = v


class PseudoPoints(BaseIO):
    """
    Holds features to be used as references / sparse points / pseudo points
    in sparse GPR methods.
    It can be populated with features computed on atomic structures using the
    extend function.

    Parameters
    ----------

    representation : Calculator
        Representation calculator associated with the kernel

    Methods
    -------
    extend(self, atoms_list, selected_indices):
        Add pseudo points to the already existing collection

        Parameters
        ----------
        atoms_list : AtomList or ManagerCollection (C++ class)
            Container of atomic structures that holds features compatible with
            representation.
    """

    def __init__(self, representation):
        self.representation = representation
        if 'SphericalInvariants' in str(representation):
            self._pseudo_points = _pseudo_points['PseudoPointsBlockSparse_SphericalInvariants'](
            )
        else:
            raise ValueError(
                'No pseudo point is appropiate for ' + str(representation))

    def get_init_params(self):
        init_params = dict(representation=self.representation)
        return init_params

    def _set_data(self, data):
        self._pseudo_points = self._pseudo_points.from_dict(
            data['pseudo_points'])

    def _get_data(self):
        return dict(pseudo_points=self._pseudo_points.to_dict())

    def extend(self, atoms_list, selected_indices):
        if isinstance(atoms_list, AtomsList):
            self._pseudo_points.extend(
                self.representation._representation, atoms_list.managers, selected_indices)
        else:
            self._pseudo_points.extend(
                self.representation._representation, atoms_list, selected_indices)

    def size(self):
        return self._pseudo_points.size()

    def get_features(self):
        return self._pseudo_points.get_features()
