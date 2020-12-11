import ase.io

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

    default_parameters = {}
    "Default parameters"

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
        self.atoms_template = ase.io.read(structure_template)
        self.manager = None

    def calculate(
        self,
        atoms=None,
        properties=["energy", "forces", "stress"],
        system_changes=all_changes,
    ):
        """ TODO: read in atomic positions and cell. Update, compute energy
        forces stress and return them in whatever format """
        Calculator.calculate(self, atoms, properties, system_changes)

        if self.manager is None:
            #  happens at the begining of the MD run
            at = self.atoms.copy()
            at.wrap(eps=1e-11)
            self.manager = [at]
        elif isinstance(self.manager, AtomsList):
            structure = unpack_ase(self.atoms, wrap_pos=True)
            structure.pop("center_atoms_mask")
            self.manager[0].update(**structure)

        self.manager = self.representation.transform(self.manager)

        energy = self.model.predict(self.manager)
        self.results["energy"] = energy
        self.results["free_energy"] = energy
        if "forces" in properties:
            self.results["forces"] = self.model.predict_forces(self.manager)
        if "stress" in properties:
            self.results["stress"] = self.model.predict_stress(self.manager).flatten()

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
