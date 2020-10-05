from ..utils import BaseIO
from ..lib import compute_sparse_kernel_gradients

import numpy as np


class KRR(BaseIO):
    """Kernel Ridge Regression model. Only supports sparse GPR
    training for the moment.

    Parameters
    ----------
    weights : np.array
        weights of the model

    kernel : Kernel
        kernel class used to train the model

    X_train : SparsePoints
        reference samples used for the training

    self_contributions : dictionary
        map atomic number to the property baseline, e.g. isolated atoms
        energies when the model has been trained on total energies.
    """

    def __init__(self, weights, kernel, X_train, self_contributions):
        # Weights of the krr model
        self.weights = weights
        self.kernel = kernel
        self.X_train = X_train
        self.self_contributions = self_contributions
        self.target_type = kernel.target_type

    def _get_property_baseline(self, managers):
        """build total baseline contribution for each prediction"""
        if self.target_type == 'Structure':
            Y0 = np.zeros(len(managers))
            for i_manager, manager in enumerate(managers):
                for center in manager:
                    Y0[i_manager] += self.self_contributions[center.atom_type]
        elif self.target_type == 'Atom':
            n_centers = 0
            for manager in managers:
                n_centers += len(manager)
            Y0 = np.zeros(n_centers)
            i_center = 0
            for manager in managers:
                for center in manager:
                    Y0[i_center] = self.self_contributions[center.atom_type]
                    i_center += 1
        return Y0

    def _preprocess_input(self, managers, compute_gradients=False):
        """compute prediction kernel and total baseline contributions"""
        kernel = self.kernel(managers, self.X_train, (compute_gradients, False))
        Y0 = self._get_property_baseline(managers)
        return kernel, Y0

    def predict(self, managers, compute_gradients=False):
        """Predict properties associated with the atomic structures in managers
        or their derivative w.r.t. atomic positions (if compute_gradients==True).

        Parameters
        ----------
        managers : AtomsList
            list of atomic structures with already computed features compatible
            with representation in kernel
        compute_gradients : bool, optional
            predict the gradients of the property w.r.t atomic positions,
            by default False

        Returns
        -------
        np.array
            predictions
        """

        if compute_gradients is False:
            KNM, Y0 = self._preprocess_input(managers, compute_gradients)
            return Y0 + np.dot(KNM, self.weights).reshape((-1))
        else:
            rep = self.kernel._representation
            gradients = compute_sparse_kernel_gradients(rep,
                        self.kernel._kernel, managers.managers,
                        self.X_train._sparse_points,
                        self.weights.reshape((1, -1)))

            return gradients


    def get_weights(self):
        return self.weights

    def _get_init_params(self):
        init_params = dict(weights=self.weights, kernel=self.kernel,
                           X_train=self.X_train,
                           self_contributions=self.self_contributions)
        return init_params

    def _set_data(self, data):
        super()._set_data(data)

    def _get_data(self):
        return super()._get_data()

    def get_representation_calculator(self):
        return self.kernel._rep


def train_gap_model(kernel, managers, KNM_, X_pseudo, y_train,
                    self_contributions, grad_train=None, lambdas=None, jitter=1e-8):
    """
    Defines the procedure to train a GAP model [1]:
    .. math::
        Y(X) = \sum_{(i,a)\in X} y_a(X_i^a),
    where :math:`Y(X)` is the predicted property function associated with the
    atomic structure :math:`X$, :math:`i` and :math:`a` are the index and
    species of the atoms in structure :math:`X` and :math:`y_a(X_i^a)` is the
    atom centered model that depends on the central atomic species.
    The individual predictions are given by:
    .. math::
        y_a(X_i) = \sum_m \alpha_m k(T_m^a, X_i^a),
    where :math:`k(\cdot,\cdot)` is a kernel function, :math:`\alpha_m` are the
    weights of the model and :math:`T_m` are the M pseudo points used to train
    the model.

    The training is given by:
    .. math::
        \bm{\alpha} =  K^{-1} \bm{Y},
    where :math:`K` is given by:
    .. math::
        K = K_{MM} + K_{MN} \Lambda^{-2} K_{NM},
    :math:`\bm{Y}=K_{MN} \Lambda^{-2} \bm{y}$, :math:`\bm{y}` the training
    targets and :math:`\Lambda` the regularization matrix.
    The regularization matrix is chossen to be diagonal:
    .. math::
        \Lambda^{-1}_{nn} = \delta_{nn} * lambdas[0] / \sigma_{\bm{y}} * np.sqrt(Natoms)
    for the targets and
    .. math::
        \Lambda^{-1}_{nn} = \delta_{nn} * lambdas[1] / \sigma_{\bm{y}},
    for the derivatives of the targets w.r.t. the atomic positions and
    :math:`\sigma_{\bm{y}}` is the standard deviation of the target property
    (not derivatives).

    Parameters
    ----------
    kernel : Kernel
        SparseKernel to compute KMM and initialize the model. It was used to
        build KNM_.
    managers : AtomsList
        Training structures with features (and gradients)
    KNM_ : np.array
        kernel matrix to use in the training, typically computed with:
        KNM = kernel(managers, X_pseudo)
        KNM_down = kernel(managers, X_pseudo, grad=(True, False))
        KNM = np.vstack([KNM, KNM_down])
        when training with derivatives.
    X_pseudo : SparsePoints
        basis samples to use in the model's interpolation
    y_train : np.array
        reference property
    self_contributions : dictionary
        map atomic number to the property baseline, e.g. training on
        total energies is not very recommended so one would provide
        the corresponding isolated atom energies so that the model
        can be trained on the corresponding formation energies and
        still predict total energies.
    grad_train : np.array, optional
        derivatives of y_train w.r.t. to the atomic motion, e.g.
        minus interatomic forces, by default None
    lambdas : list/tuple, optional
        regularisation parameter for the training, i.e. lambdas[0] -> property
        and lambdas[1] -> gradients of the property, by default None
    jitter : double, optional
        small jitter for the numerical stability of solving the linear system,
        by default 1e-8

    Returns
    -------
    KRR
        a trained model that can predict the property and its gradients

    .. [1] Ceriotti, M., Willatt, M. J., & Csányi, G. (2018).
        Machine Learning of Atomic-Scale Properties Based on Physical Principles.
        In Handbook of Materials Modeling (pp. 1–27). Springer, Cham.
        https://doi.org/10.1007/978-3-319-42913-7_68-1
    """
    KMM = kernel(X_pseudo)
    Y = y_train.reshape((-1, 1)).copy()
    KNM = KNM_.copy()
    n_centers = Y.shape[0]
    Natoms = np.zeros(n_centers)
    Y0 = np.zeros((n_centers, 1))
    for iframe, manager in enumerate(managers):
        Natoms[iframe] = len(manager)
        for center in manager:
            Y0[iframe] += self_contributions[center.atom_type]
    Y = Y - Y0
    delta = np.std(Y)
    # lambdas[0] is provided per atom hence the '* np.sqrt(Natoms)'
    # the first n_centers rows of KNM are expected to refer to the
    # property
    KNM[:n_centers] /= lambdas[0] / delta * np.sqrt(Natoms)[:, None]
    Y /= lambdas[0] / delta * np.sqrt(Natoms)[:, None]

    if grad_train is not None:
        KNM[n_centers:] /= lambdas[1] / delta
        F = grad_train.reshape((-1, 1)).copy()
        F /= lambdas[1] / delta
        Y = np.vstack([Y, F])

    KMM[np.diag_indices_from(KMM)] += jitter

    K = KMM + np.dot(KNM.T, KNM)
    Y = np.dot(KNM.T, Y)
    weights = np.linalg.lstsq(K, Y, rcond=None)[0]
    model = KRR(weights, kernel, X_pseudo, self_contributions)

    # avoid memory clogging
    del K, KMM
    K, KMM = [], []

    return model
