from ..utils import BaseIO
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

    def __init__(self, model, representation, **kwargs):
        """ TODO: get in the name of the model, and do all of the initialization stuff that is needed. 
        If atomic info is needed, it will have to be read in here, e.g. given a reference xyz file.  """
        super(ASEMLCalculator, self).__init__(**kwargs)
        self.model = model
        self.representation = representation
        self.kwargs = kwargs
        self.manager = None

    def calculate(
        self,
        atoms=None,
        properties=["energy", "forces", "stress"],
        system_changes=all_changes,
    ):
        """ TODO: read in atomic positions and cell. Update, compute energy forces stress and return them in whatever format """
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
        init_params = dict(model=self.model, representation=self.representation)
        init_params.update(**self.kwargs)
        return init_params

    def _set_data(self, data):
        self.manager = None
        super()._set_data(data)

    def _get_data(self):
        return super()._get_data()
