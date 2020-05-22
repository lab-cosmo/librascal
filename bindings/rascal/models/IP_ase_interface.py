from ..utils import BaseIO
from ase.calculators.calculator import Calculator, all_changes
from copy import deepcopy


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

    implemented_properties = ["energy", "forces"]
    "Properties calculator can handle (energy, forces, ...)"

    default_parameters = {}
    "Default parameters"

    nolabel = True

    def __init__(self, model, representation, **kwargs):
        super(ASEMLCalculator, self).__init__(**kwargs)
        self.model = model
        self.representation = representation
        self.kwargs = kwargs

    def calculate(
        self, atoms=None, properties=["energy"], system_changes=all_changes
    ):
        Calculator.calculate(self, atoms, properties, system_changes)

        X = [self.atoms]
        managers = self.representation.transform(X)

        energy = self.model.predict(managers)
        forces = -self.model.predict(managers, compute_gradients=True)

        self.results["energy"] = energy
        self.results["free_energy"] = energy
        self.results["forces"] = forces

    def _get_init_params(self):
        init_params = dict(
            model=self.model, representation=self.representation
        )
        init_params.update(**self.kwargs)
        return init_params

    def _set_data(self, data):
        super()._set_data(data)

    def _get_data(self):
        return super()._get_data()
