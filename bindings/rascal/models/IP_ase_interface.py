from ase.calculators.calculator import Calculator, all_changes
from copy import deepcopy


class ASEMLCalculator(Calculator):
    """Wrapper class to use a rascal model as an interatomic potential in ASE

    Parameters
    ----------
    model : class
        a trained model of the rascal library that can predict the energy and
        derivaties of the energy w.r.t. atomic positions
    representation : class
        a representation calculator of rascal compatible with the trained model
    """

    implemented_properties = ['energy', 'forces']
    'Properties calculator can handle (energy, forces, ...)'

    default_parameters = {}
    'Default parameters'

    nolabel = True

    def __init__(self, model, representation, **kwargs):
        Calculator.__init__(self, **kwargs)

        self.model = model
        self.representation = representation

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        X = [self.atoms]
        managers = self.representation.transform(X)

        energy = self.model.predict(managers)
        forces = -self.model.predict(managers, compute_gradients=True)

        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces
