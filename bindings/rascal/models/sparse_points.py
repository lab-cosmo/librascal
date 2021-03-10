from ..utils import BaseIO
from ..lib._rascal.models import kernels
from ..neighbourlist import AtomsList

# names of existing pseudo points implementation on the pybinding side.
_sparse_points = {}
for k, v in kernels.__dict__.items():
    if "SparsePoints" in k:
        name = k
        _sparse_points[name] = v


class SparsePoints(BaseIO):
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
        Extends the list of sparse points with the features contained in
        atoms_list (AtomList class a wrapper for the ManagerCollection class in
        c++) using selected_indices to choose which atom centered environment
        will be add to the sparse points.

        Parameters
        ----------
        atoms_list : AtomList or ManagerCollection (C++ class)
            Container of atomic structures that holds features compatible with
            representation.

        selected_indices : a list (same length as atoms_list) containing the
            list of selected atomic environment indices within their respective
            structures, e.g. [[0,3,7],[]] corresponds to selecting the
            representation of environment 0, 3 and 7 in the first structure of
            atoms_list and none in the second structure.
    """

    def __init__(self, representation):
        super(SparsePoints, self).__init__()
        self.representation = representation
        if "SphericalInvariants" in str(representation):
            self._sparse_points = _sparse_points[
                "SparsePointsBlockSparse_SphericalInvariants"
            ]()
        else:
            raise ValueError("No pseudo point is appropiate for " + str(representation))

    def _get_init_params(self):
        init_params = dict(representation=self.representation)
        return init_params

    def _set_data(self, data):
        super()._set_data(data)
        self._sparse_points = self._sparse_points.from_dict(data["sparse_points"])

    def _get_data(self):
        data = super()._get_data()
        data.update(sparse_points=self._sparse_points.to_dict())
        return data

    def extend(self, atoms_list, selected_indices):
        if isinstance(atoms_list, AtomsList):
            self._sparse_points.extend(
                self.representation._representation,
                atoms_list.managers,
                selected_indices,
            )
        else:
            self._sparse_points.extend(
                self.representation._representation,
                atoms_list,
                selected_indices,
            )

    def size(self):
        return self._sparse_points.size()

    def __len__(self):
        return self.size()

    def get_features(self):
        return self._sparse_points.get_features()
