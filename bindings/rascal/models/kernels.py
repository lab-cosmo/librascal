from ..lib._rascal.models.kernels import Kernel as Kernelcpp
from ..neighbourlist import AtomsList
import json


class Kernel(object):

    """ Computes the kernel from an representation """

    def __init__(self, representation, name='Cosine', target_type='structure',
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
            Type of target (prediction) properties, either 'atomic' (kernel is
            between atomic environments) or 'structure' (kernel is summed over
            atoms in a structure) (default)

        """
        hypers = dict(name=name,target_type=target_type)
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
