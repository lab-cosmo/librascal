from ..lib._rascal.models.kernels import Kernel as Kernelcpp
from ..neighbourlist import AtomsList
import json


class Kernel(object):

    """Computes the kernel from an representation
    Initialize the kernel with the given parameters

    Parameters
    ----------

    representation : Calculator

        Representation calculator associated with the kernel

    name : string

        Type of kernel, 'Cosine' (aka dot-product) is the default and
        (currently) only option.

    target_type : string

        Type of target (prediction) properties, must be either 'Atom' (the
        kernel is between atomic environments) or 'Structure' (the kernel is
        summed over atoms in a structure), which is the default

    Methods
    -------
    __call__(X, Y=None)
        Compute the kernel.

        Parameters
        ----------
        X : AtomList or ManagerCollection (C++ class)

            Container of atomic structures.

        Returns
        -------
        kernel_matrix: ndarray

    """

    def __init__(self, representation, name='Cosine', target_type='Structure',
                 **kwargs):

        # This case cannot be handled by the C++ side because the type cannot be
        # deduced from arguments inside a json. It therefore has to be cast on
        # the C++ side and \zeta has to be checked here.
        if (name == 'Cosine' and 'zeta' in kwargs):
            # should be positive integer
            zeta = kwargs['zeta']

            if not(zeta > 0 and zeta % 1 == 0):
                raise ValueError(
                    "The given zeta has to be a positive integer.")

        hypers = dict(name=name, target_type=target_type)
        hypers.update(**kwargs)
        hypers_str = json.dumps(hypers)
        self._representation = representation._representation

        self._kernel = Kernelcpp(hypers_str)

    def __call__(self, X, Y=None):
        if isinstance(X, AtomsList):
            X = X.managers

        if Y is None:
            return self._kernel.compute(self._representation, X)

        else:
            if isinstance(Y, AtomsList):
                Y = Y.managers

            return self._kernel.compute(self._representation, X, Y)
