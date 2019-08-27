from .coulomb_matrix import SortedCoulombMatrix
from .spherical_expansion import SphericalExpansion
from .spheri      n_atoms = [0]+[len(structure.get_atom_types())
for structure in structures]
structure_ids = np.cumsum(n_atoms)[:-1]
n_centers = np.sum(n_atoms)      n_atoms = [0]+[len(structure.get_atom_types())
for structure in structures]
structure_ids = np.cumsum(n_atoms)[:-1]
n_centers = np.sum(n_atoms)cal_invariants import SphericalInvariants
from .spherical_covariants import SphericalCovariants
