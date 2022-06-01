from ..utils import BaseIO
from ..lib._rascal.models import (
    compute_numerical_kernel_gradients as _compute_numerical_kernel_gradients,
)
from ..lib._rascal.models.kernels import Kernel as Kernelcpp
from ..lib._rascal.models.kernels import SparseKernel as SparseKernelcpp
from ..neighbourlist import AtomsList
from ..utils import BaseIO
from .sparse_points import SparsePoints
import json


import numpy as np


class Kernel(BaseIO):

    """
    Initialize the kernel with the given representation and parameters

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
        or 'Sparse' (computing GAP [2]_ like kernel for sparse kernel methods like
        Subset of Regressors)

    Notes
    -----
    In the following
    we refer to the training samples with 'N' and, in the case of
    sparse kernels [1]_, we refer to the pseudo points with 'M'. So a
    kernel between the training samples and the pseudo points is
    'KNM'. For more information on sparse kernels see
    :meth:`rascal.models.krr.train_gap_model`.

    .. [1] Joaquin Quiñonero-Candela, Carl Edward Rasmussen;
            A Unifying View of Sparse Approximate Gaussian Process Regression,
            6(Dec):1939--1959, 2005.
            http://www.jmlr.org/papers/v6/quinonero-candela05a.html

    .. [2] Ceriotti, M., Willatt, M. J., & Csányi, G. (2018).
    Machine Learning of Atomic-Scale Properties Based on Physical Principles.
    In Handbook of Materials Modeling (pp. 1–27). Springer, Cham.
    https://doi.org/10.1007/978-3-319-42913-7_68-1

    """

    def __init__(
        self,
        representation,
        name="Cosine",
        kernel_type="Full",
        target_type="Structure",
        **kwargs,
    ):
        # This case cannot be handled by the c++ side because c++ cannot deduce the
        # type from arguments inside a json, so it has to be casted in the c++
        # side. Therefore zeta has to be checked here.
        if name == "Cosine" and "zeta" in kwargs:
            # should be positive integer
            zeta = kwargs["zeta"]
            if not (zeta > 0 and isinstance(zeta, int)):
                raise ValueError("The given zeta has to be a positive integer.")
        elif name == "GAP" and "zeta" in kwargs:
            # should be positive integer
            zeta = kwargs["zeta"]
            if not (zeta > 0 and isinstance(zeta, int)):
                raise ValueError("The given zeta has to be a positive integer.")
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
        if "Sparse" in kernel_type:
            self._kernel = SparseKernelcpp(hypers)
        else:
            self._kernel = Kernelcpp(hypers)

    def _get_init_params(self):
        init_params = dict(
            representation=self._rep,
            name=self.name,
            kernel_type=self.kernel_type,
            target_type=self.target_type,
        )
        init_params.update(**self._kwargs)
        return init_params

    def _set_data(self, data):
        super()._set_data(data)

    def _get_data(self):
        return super()._get_data()

    def __call__(self, X, Y=None, grad=(False, False), compute_neg_stress=False):
        """
        Compute the kernel.

        Parameters
        ----------
            X : AtomList or ManagerCollection (C++ class)
                Container of atomic structures.

            Y : AtomList, ManagerCollection or SparsePoints* (C++ class).

            grad : tuple specifying if the kernel should be computed using
                gradient
                of the representation w.r.t. the atomic positions, e.g. (True,
                False) corresponds to the gradient of the 1st argument of the
                kernel w.r.t. the atomic positions.

            compute_neg_stress : if gradients are computed and True then compute
                also the kernel associated with the stress in Voigt format.

        Returns
        -------
            kernel_matrix: ndarray
        """
        if isinstance(X, AtomsList):
            X = X.managers
        if Y is None and grad == (False, False):
            # compute a kernel between features and themselves
            if self.kernel_type == "Full":
                return self._kernel.compute(self._representation, X)
            elif self.kernel_type == "Sparse":
                if isinstance(X, SparsePoints):
                    X = X._sparse_points
                # compute the KMM matrix
                return self._kernel.compute(X)
        elif grad == (True, False) and "Sparse" in self.kernel_type:
            if isinstance(Y, SparsePoints):
                Y = Y._sparse_points
            # compute the block of the KNM matrix corresponding to forces
            return self._kernel.compute_derivative(
                self._representation, X, Y, compute_neg_stress
            )
        elif grad == (False, False):
            # compute the kernel between two sets of features
            if isinstance(Y, AtomsList):
                # to make predictions with the full covariance
                Y = Y.managers
            elif isinstance(Y, SparsePoints):
                # to make predictions with a sparse kernel method
                Y = Y._sparse_points
            return self._kernel.compute(self._representation, X, Y)
        else:
            raise NotImplementedError(
                "The configuration: {} is not implemented for kernel {} in {} mode.".format(
                    grad, self.name, self.kernel_type
                )
            )


class KernelDirect(BaseIO):

    """Compute a kernel directly by specifying the feature matrices

    This class aims to be compatible with the KRR class and gaptools
    infrastructure, while not requiring a link to a librascal representation
    calculator.  The tradeoff is lower computational efficiency, especially
    for sparse multi-species systems.

    TODO implement structure and gradient kernels
    """

    def __init__(
        self,
        kernel_power,
        representation=None,
        kernel_name="Cosine",
        target_type="Atom",
        features_key="features",
    ):
        """Make a kernel function with the given parameters

        Currently only supports atom-wise dot-product (aka Cosine)
        kernels, raised to some power.

        Parameters
        ----------
        kernel_power : int or float
            Exponent of the polynomial kernel, sometimes called zeta

        Other parameters
        ----------------
        representation : object
            Optional representation calculator, currently unused
        kernel_name : str
            Type of kernel; must be "Cosine"
        target_type : str
            Selects whether the kernel is computed for atoms ('Atom')
            or structures ('Structure').
        features_key : str
            For structure kernels, which key is used to access the
            features in the provided ASE Atoms objects.
            This can be modified in the kernel __call__() as well.

        Returns
        -------
        kernel : function
            Function for computing the kernel between two feature matrices
        """
        if kernel_name != "Cosine":
            raise ValueError("Only cosine kernels are supported")
        valid_target_types = ["Atom", "Structure"]
        if target_type not in valid_target_types:
            raise ValueError(f"'target_type' must be one of: {valid_target_types!s}")
        self.kernel_name = kernel_name
        self.target_type = target_type
        self.features_key = features_key
        self.kernel_power = kernel_power
        self._representation = representation

    def _self_kernel(self, features):
        """Compute the kernel of the feature matrix with itself"""
        return (features @ features.T) ** self.kernel_power

    def __call__(self, features, other_features=None, grad=False, features_key=None):
        """Compute the kernel between the feature matrices

        Parameters
        ----------
        features : 2-D array or list(ase.Atoms)
            First set of features.  If target_type=='Structure', should
            instead be a list of ASE Atoms with the features for each
            structure stored in the Atoms.arrays attribute, e.g. using
            librascal.neighbourlist.store_features_ase_atoms()
        other_features : 2-D array, optional
            Second set of features (e.g. sparse points)
            If not provided, the kernel is computed between the
            first feature set and itself.
        features_key : str
            For structure kernels, the dictionary key used to access
            features in the ASE Atoms objects.  If provided, overrides
            the key provided in the constructor.
        grad : bool
            Whether to compute the _gradient_ of the kernel with respect
            to atomic positions.  This operation is only applied to the
            first set of features, since the double-gradient kernel is
            practically never used.  Not supported for the self-kernel
            (i.e. if other_features is None).

        Returns
        -------
        kernel : 2-D array
            The kernel between the requested features
        """
        # Compatibility with the original Kernel call function
        if isinstance(grad, tuple):
            grad = grad[0]
        if other_features is None:
            if grad:
                raise ValueError("Gradients are not supported for the self-kernel")
            return self._self_kernel(features)
        else:
            if grad:
                # TODO implement
                raise NotImplemented("Gradient kernel WIP, sorry")
            else:
                if self.target_type == "Structure":
                    output_kernel = np.empty((len(features), other_features.shape[0]))
                    # This may be a somewhat slow way of doing the computation,
                    # but avoids hogging large amounts of memory at once.
                    # Consider implementing a "blocking" option to compromise the
                    # two needs in the future.
                    offset = 0
                    if features_key is None:
                        features_key = self.features_key
                    for structure_idx, structure in enumerate(features):
                        try:
                            structure_features = structure.arrays[features_key]
                        except AttributeError as ae:
                            raise ValueError(
                                "Need a list of ASE Atoms for structure kernels"
                            ) from ae
                        except KeyError as ke:
                            raise ValueError(
                                f"No features found under key '{features_key}'"
                            ) from ke
                        output_kernel[structure_idx] = np.sum(
                            (structure_features @ other_features.T)
                            ** self.kernel_power,
                            axis=0,
                        )
                    return output_kernel
                else:
                    return (features @ other_features.T) ** self.kernel_power

    def _set_data(self, data):
        super()._set_data(data)

    def _get_data(self):
        return super()._get_data()

    def _get_init_params(self):
        return dict(
            kernel_power=self.kernel_power,
            kernel_name=self.kernel_name,
            target_type=self.target_type,
            representation=self._representation,
        )


def compute_numerical_kernel_gradients(
    kernel, calculator, managers, sparse_points, h_disp, compute_neg_stress=True
):
    """This function is used for testing the numerical kernel gradient"""
    return _compute_numerical_kernel_gradients(
        kernel._kernel,
        calculator._representation,
        managers.managers,
        sparse_points._sparse_points,
        h_disp,
        compute_neg_stress,
    )
