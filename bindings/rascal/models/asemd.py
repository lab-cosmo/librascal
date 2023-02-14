from ..utils import BaseIO
from ase.calculators.calculator import Calculator, all_changes
from copy import deepcopy
from ..neighbourlist.structure_manager import AtomsList, unpack_ase


class ASEMLCalculator(Calculator, BaseIO):
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
        super(ASEMLCalculator, self).__init__(**kwargs)
        if not isinstance(model, list):
            model = [model]

        if not isinstance(representation,list):
            representation = [representation]

        if len(model)!=len(representation):
            raise ValueError('should provide a representation per model. len(model)!=len(representation)')

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
        Calculator.calculate(self, atoms, properties, system_changes)

        if self.manager is None:
            #  happens at the begining of the MD run
            at = self.atoms.copy()
            at.wrap(eps=1e-11)
            self.manager = [at]*len(self.model)
        elif isinstance(self.manager, list) and isinstance(self.manager[0], AtomsList):
            structure = unpack_ase(self.atoms, wrap_pos=True)
            structure.pop("center_atoms_mask")
            for i in range(len(self.manager)):
                self.manager[i][0].update(**structure)

        energy = 0
        for i in range(len(self.manager)):
            self.manager[i] = self.representation[i].transform(self.manager[i])
            energy += self.model[i].predict(self.manager[i])
            if "forces" in properties:
                if i == 0:
                    self.results["forces"] = self.model[i].predict_forces(self.manager[i])
                else:
                    self.results["forces"] += self.model[i].predict_forces(self.manager[i])
            if "stress" in properties:
                if i == 0:
                    self.results["stress"] = self.model[i].predict_stress(self.manager[i]).flatten()
                else:
                    self.results["stress"] += self.model[i].predict_stress(self.manager[i]).flatten()

        self.results["energy"] = energy
        self.results["free_energy"] = energy

    def _get_init_params(self):
        init_params = dict(model=self.model[0], representation=self.representation[0])
        init_params.update(**self.kwargs)
        return init_params

    def _set_data(self, data):
        self.manager = None
        super()._set_data(data)

    def _get_data(self):
        return super()._get_data()
