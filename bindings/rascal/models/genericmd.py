import ase.io
import ase.units
import numpy as np

from ..utils import BaseIO, load_obj
from ..neighbourlist.structure_manager import AtomsList, unpack_ase


class GenericMDCalculator(BaseIO):

    """Generic MD driver for a librascal model

    Initialize with model JSON and a structure template, and calculate
    energies and forces based on position/cell updates _assuming the
    order and identity of atoms does not change_.
    """

    def __init__(self, model_json, structure_template, assume_pbc=True):
        """Initialize a model and structure template

        Parameters
        ----------
        model_json  Filename for a JSON file defining the potential
        structure_template
                    Filename for an ASE-compatible Atoms object, used
                    only to initialize atom types and numbers
        assume_pbc  Force the PBC to True (to avoid subtle issues that
                    sometimes arise when no cell is defined), default True
        """
        super(GenericMDCalculator, self).__init__()
        self.model_filename = model_json
        self.model = load_obj(model_json)
        self.representation = self.model.get_representation_calculator()
        self.template_filename = structure_template
        self.atoms = ase.io.read(structure_template, 0)
        if assume_pbc:
            self.atoms.pbc = True
        self.manager = None
        self.matrix_indices_in_voigt_notation = [
            (0, 0),
            (1, 1),
            (2, 2),
            (1, 2),
            (0, 2),
            (0, 1),
        ]

    def calculate(self, positions, cell_matrix):
        """Calculate energies and forces from position/cell update

        positions   Atomic positions (Nx3 matrix)
        cell_matrix Unit cell (in ASE format, cell vectors as rows)

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

    def _get_init_params(self):
        init_params = dict(
            model_json=self.model_filename, structure_template=self.template_filename
        )
        return init_params

    def _set_data(self, data):
        self.manager = None
        super()._set_data(data)

    def _get_data(self):
        return super()._get_data()
