from ..lib._rascal.models.kernels import Kernel as Kernelcpp
from ..lib._rascal.models.kernels import SparseKernel as SparseKernelcpp
from ..neighbourlist import AtomsList
from .pseudo_points import PseudoPoints
from ..utils import return_deepcopy, BaseIO
import json


class Kernel(BaseIO):

    """
    Computes the kernel from an representation
    Initialize the kernel with the given parameters

    Parameters
    ----------

    representation : Calculator
        Representation calculator associated with the kernel

    name : string
        Type of kernel, 'Cosine' (aka dot-product) is the default and
        (currently) only option.

    target_type : string
        Type of target (prediction) properties, must be either 'Atom' (the kernel
        is between atomic environments) or 'Structure' (the kernel is summed over
        atoms in a structure), which is the default

    kernel_type : string
        Type of kernel method, either 'Full' (computing exact covariance matrix)
        or 'Sparse' (computing GAP like kernel for sparse kernel methods like
        Subset of Regressors)

    Methods
    -------
    __call__(X, Y=None, grad=(False, False))
        Compute the kernel.

        Parameters
        ----------
        X : AtomList or ManagerCollection (C++ class)
            Container of atomic structures.

        Y : AtomList, ManagerCollection or PseudoPoints* (C++ class).

        grad : tuple specifying if the kernel should be computed using gradient
               of the representation w.r.t the atomic positions.

        Returns
        -------
        kernel_matrix: ndarray

    """

    def __init__(self, representation, name='Cosine', kernel_type='Full', target_type='Structure',
                 **kwargs):

        # This case cannot handled by the c++ side because c++ cannot deduce the
        # type from arguments inside a json, so it has to be casted in the c++
        # side. Therefore zeta has to be checked here.
        if (name == 'Cosine' and 'zeta' in kwargs):
            # should be positive integer
            zeta = kwargs['zeta']
            if not(zeta > 0 and zeta % 1 == 0):
                raise ValueError(
                    "The given zeta has to be a positive integer.")
        elif (name == 'GAP' and 'zeta' in kwargs):
            # should be positive integer
            zeta = kwargs['zeta']
            if not(zeta > 0 and zeta % 1 == 0):
                raise ValueError(
                    "The given zeta has to be a positive integer.")
        else:
            raise RuntimeError("Kernel name must be one of: Cosine, GAP.")
        hypers = dict(name=name, target_type=target_type)
        hypers.update(**kwargs)
        self._rep = representation
        self._representation = representation._representation
        self.name = name
        self._kwargs = kwargs
        self.kernel_type = kernel_type
        self.target_type = target_type
        if 'Sparse' in kernel_type:
            self._kernel = SparseKernelcpp(hypers)
        else:
            self._kernel = Kernelcpp(hypers)

    def get_init_params(self):
        init_params = dict(representation=self._rep,
                           name=self.name,
                           kernel_type=self.kernel_type,
                           target_type=self.target_type)
        init_params.update(**self._kwargs)
        return init_params

    def _set_data(self, data):
        pass

    def _get_data(self):
        return dict()

    def __call__(self, X, Y=None, grad=(False, False)):
        if isinstance(X, AtomsList):
            X = X.managers
        if Y is None and grad == (False, False):
            # compute a kernel between features and themselves
            if self.kernel_type == 'Full':
                return self._kernel.compute(self._representation, X)
            elif self.kernel_type == 'Sparse':
                if isinstance(X, PseudoPoints):
                    X = X._pseudo_points
                # compute the KMM matrix
                return self._kernel.compute(X)
        elif grad == (True, False) and 'Sparse' in self.kernel_type:
            if isinstance(Y, PseudoPoints):
                Y = Y._pseudo_points
            # compute the block of the KNM matrix corresponding to forces
            return self._kernel.compute_derivative(self._representation, X, Y)
        elif grad == (False, False):
            # compute the kernel between two sets of features
            if isinstance(Y, AtomsList):
                # to make predictions with the full covariance
                Y = Y.managers
            elif isinstance(Y, PseudoPoints):
                # to make predictions with a sparse kernel method
                Y = Y._pseudo_points
            return self._kernel.compute(self._representation, X, Y)
        else:
            raise NotImplementedError('The configuration: {} is not implemented for kernel {} in {} mode.'.format(
                grad, self.name, self.kernel_type))
