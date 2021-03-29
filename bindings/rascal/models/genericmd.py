"""Generic calculator-style interface for MD"""
import ase.io
import numpy as np

from ..utils import BaseIO, load_obj
from ..neighbourlist.structure_manager import AtomsList, unpack_ase


class GenericMDCalculator:

    """Generic MD driver for a librascal model

    Initialize with model JSON and a structure template, and calculate
    energies and forces based on position/cell updates _assuming the
    order and identity of atoms does not change_.
    """

    matrix_indices_in_voigt_notation = (
        (0, 0),
        (1, 1),
        (2, 2),
        (1, 2),
        (0, 2),
        (0, 1),
    )

    def __init__(
        self, model_json, is_periodic, structure_template=None, atomic_numbers=None
    ):
        """Initialize a model and structure template

        Parameters
        ----------
        model_json  Filename for a JSON file defining the potential
        is_periodic Specify whether the simulation is periodic or not
                    This helps avoid confusion if a geometry's "periodic"
                    flags have been set improperly, which can happen e.g.
                    if ASE cannot read the cell information in a file.  If
                    using a structure template and this is set to True,
                    will raise an error unless at least one of the PBC
                    flags in the structure template is on.  If set to
                    False, will raise an error if all PBC flags are not
                    off.  Set to None to skip PBC checking.  If not using a
                    structure template, this setting will determine the PBC
                    flags of the created atomic structure.
        structure_template
                    Filename for an ASE-compatible Atoms object, used
                    only to initialize atom types and numbers
        atomic_numbers
                    List of atom types (atomic numbers) to initialize
                    the atomic structure in case no structure template
                    is given
        """
        super(GenericMDCalculator, self).__init__()
        self.model_filename = model_json
        self.model = load_obj(model_json)
        self.representation = self.model.get_representation_calculator()
        self.manager = None
        # Structure initialization
        self.is_periodic = is_periodic
        if structure_template is not None:
            self.template_filename = structure_template
            self.atoms = ase.io.read(structure_template, 0)
            if (is_periodic is not None) and (
                is_periodic != np.any(self.atoms.get_pbc())
            ):
                raise ValueError(
                    "Structure template PBC flags: "
                    + str(self.atoms.get_pbc())
                    + " incompatible with 'is_periodic' setting"
                )
        elif atomic_numbers is not None:
            self.atoms = ase.Atoms(numbers=atomic_numbers, pbc=is_periodic)
        else:
            raise ValueError(
                "Must specify one of 'structure_template' or 'atomic_numbers'"
            )

    def calculate(self, positions, cell_matrix):
        """Calculate energies and forces from position/cell update

        positions   Atomic positions (Nx3 matrix)
        cell_matrix Unit cell (in ASE format, cell vectors as rows)
                    (set to zero for non-periodic simulations)

        The units of positions and cell are determined by the model JSON
        file; for now, only Å is supported.  Energies, forces, and
        stresses are returned in the same units (eV and Å supported).

        Returns a tuple of energy, forces, and stress - forces are
        returned as an Nx3 array and stresses are returned as a 3x3 array

        Stress convention: The stresses have units eV/Å^3
        (volume-normalized) and are defined as the gradients of the
        energy with respect to the cell parameters.
        """
        # Update ASE Atoms object (we only use ASE to handle any
        # re-wrapping of the atoms that needs to take place)
        self.atoms.set_cell(cell_matrix)
        self.atoms.set_positions(positions)

        # Quick consistency checks
        if positions.shape != (len(self.atoms), 3):
            raise ValueError(
                "Improper shape of positions (is the number of atoms consistent?)"
            )
        if cell_matrix.shape != (3, 3):
            raise ValueError("Improper shape of cell info (expected 3x3 matrix)")

        # Convert from ASE to librascal
        if self.manager is None:
            #  happens at the begining of the MD run
            at = self.atoms.copy()
            at.wrap(eps=1e-11)
            self.manager = [at]
        elif isinstance(self.manager, AtomsList):
            structure = unpack_ase(self.atoms, wrap_pos=True)
            structure.pop("center_atoms_mask")
            self.manager[0].update(**structure)

        # Compute representations and evaluate model
        self.manager = self.representation.transform(self.manager)
        energy = self.model.predict(self.manager)
        forces = self.model.predict_forces(self.manager)
        stress_voigt = self.model.predict_stress(self.manager)
        stress_matrix = np.zeros((3, 3))
        stress_matrix[tuple(zip(*self.matrix_indices_in_voigt_notation))] = stress_voigt
        # Symmetrize the stress matrix (replicate upper-diagonal entries)
        stress_matrix += np.triu(stress_matrix).T
        return energy, forces, stress_matrix
