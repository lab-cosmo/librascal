from ..lib._rascal.models.kernels import Kernel as Kernelcpp
from ..neighbourlist import AtomsList
import json


class Kernel(object):

    """ Computes the kernel from an representation """

    def __init__(self, representation, name='Cosine', target_type='Structure',
                 **kwargs):
        """
        Initialize the kernel with the given parameters

        Parameters
        ----------

        representation : Calculator
            Representation calculator associated with the kernel

        name : string
            Type of kernel, 'Cosine' (aka dot-product) is the default and
            (currently) only option

        target_type : string
            Type of target (prediction) properties, either 'Atom' (kernel is
            between atomic environments) or 'Structure' (kernel is summed over
            atoms in a structure) (default)

        """

        # This case cannot handled by the c++ side because c++ cannot deduce the
        # type from arguments inside a json, so it has to be casted in the c++
        # side. Therefore zeta has to be checked here.
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
        """
        Compute the kernel.

        Parameters
        ----------
        X : AtomList or ManagerCollection (C++ class)
            Container of atomic structures.

        Returns
        -------
        kernel_matrix: ndarray
        """
        if Y is None:
            Y = X
        if isinstance(X, AtomsList):
            X = X.managers
        if isinstance(Y, AtomsList):
            Y = Y.managers
        return self._kernel.compute(self._representation, X, Y)
