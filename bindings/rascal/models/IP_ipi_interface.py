import ase.io
import ase.units

from ..utils import BaseIO, load_obj
from copy import deepcopy
from ..neighbourlist.structure_manager import AtomsList, unpack_ase


class IPICalculator(BaseIO):
    """Wrapper class to use a rascal model as an interatomic potential in ASE

    Parameters
    ----------
    model : class
        a trained model of the rascal library that can predict the energy and
        derivaties of the energy w.r.t. atomic positions
    representation : class
        a representation calculator of rascal compatible with the trained model
    """

    implemented_properties = ["energy", "forces", "stress"]
    "Properties calculator can handle (energy, forces, ...)"

    # I don't think this is needed...
    default_parameters = {}
    "Default parameters"

    # Is this something required by ASE...? If so, take it out.
    nolabel = True

    def __init__(self, model_json, structure_template):
        """Initialize a model from i-PI

        Parameters:
            model_json  Filename for a JSON file defining the potential
            structure_template
                        Filename for an ASE-compatible Atoms object, used
                        only to initialize atom types and numbers
        """
        super(IPICalculator, self).__init__()
        self.model_filename = model_json
        self.model = load_obj(model_json)
        self.representation = self.model.get_representation_calculator()
        self.template_filename = structure_template
        self.atoms = ase.io.read(structure_template)
        self.manager = None
        self.matrix_indices_in_voigt_notation = [
            (0, 0),
            (1, 1),
            (2, 2),
            (1, 2),
            (0, 2),
            (0, 1),
        ]

    def calculate(self, positions, cell_tuple):
        """Calculate energies and forces from i-PI update

        positions   Atomic positions, in atomic units (Bohr)
        cell_tuple  Unit cell (and inverse cell), in atomic units

        WARNING positions and cell are passed in atomic units.  ASE uses
        eV and Å, so we're using that in librascal by default.
        """

        # Update ASE Atoms object (probably unnecessary to use ASE as
        # intermediary, maybe we can just use the internal format...)

        # The ASE cell convention is transposed from the i-PI convention,
        # which adheres to the usual "cell vectors are columns" convention
        cell_ase = cell_tuple[0].T * ase.units.Bohr
        self.atoms.set_cell(cell_ase)
        positions_ase = positions.reshape((-1, 3)) * ase.units.Bohr
        self.atoms.set_positions(positions_ase)

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
        energy = self.model.predict(self.manager) / ase.units.Hartree
        forces = (self.model.predict_forces(self.manager).flatten()
                  / ase.units.Hartree * ase.units.Bohr)
        # TODO figure out virial sign convention
        stress_voigt = self.model.predict_stress(self.manager)
        virial_voigt_ipi = stress_voigt * self.atoms.get_volume() / ase.units.Hartree
        virial_ipi = np.zeros((3, 3))
        virial_ipi[self.matrix_indices_in_voigt_notation] = virial_voigt_ipi
        return energy, forces, virial_ipi

    def _get_init_params(self):
        init_params = dict(
            model_json=self.model_filename,
            structure_template=self.template_filename
        )
        return init_params

    def _set_data(self, data):
        self.manager = None
        super()._set_data(data)

    def _get_data(self):
        return super()._get_data()
