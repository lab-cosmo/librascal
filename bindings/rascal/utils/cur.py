from ..models import PseudoPoints
from ..utils import BaseIO

import numpy as np
from scipy.sparse.linalg import svds


def do_CUR(X, Nsel, act_on='sample', is_deterministic=False, seed=10, verbose=True):
    """ Apply CUR selection [1] of Nsel rows or columns of the
    given feature matrix X[n_samples, n_features].

    .. [1] Mahoney, M. W., & Drineas, P. (2009). CUR matrix decompositions for improved data analysis. Proceedings of the National Academy of Sciences,106(3), 697–702. https://doi.org/10.1073/pnas.0803205106
    """
    U, _, VT = svds(X, Nsel)
    if 'sample' in act_on:
        weights = np.mean(np.square(U), axis=1)
    elif 'feature' in act_on:
        weights = np.mean(np.square(VT), axis=0)
    if is_deterministic is True:
        # sorting is smallest to largest hence the minus
        sel = np.argsort(-weights)[:Nsel]
    elif is_deterministic is False:
        np.random.seed(seed)
        # sorting is smallest to largest hence the minus
        sel = np.argsort(np.random.rand(*weights.shape) - weights)[:Nsel]

    if verbose is True:
        if 'sample' in act_on:
            C = X[sel, :]
        elif 'feature' in act_on:
            C = X[:, sel]
        Cp = np.linalg.pinv(C)
        err = np.sqrt(np.sum((X - np.dot(np.dot(X, Cp), C))**2))
        print('Reconstruction RMSE={:.3e}'.format(err))

    return sel


class CURFilter(BaseIO):
    """CUR decomposition to select samples or features in a given feature matrix. Wrapper around the do_CUR function for convenience.

    Parameters
    ----------
    representation : Calculator
        Representation calculator associated with the kernel

    Nselect: int
        number of points to select. if act_on='sample per specie' then it should
        be a dictionary mapping atom type to the number of samples, e.g.
        Nselect = {1:200,6:100,8:50}.

    act_on: string
        Select how to apply the selection. Can be either of 'sample',
        'sample per species','feature'.
        For the moment only 'sample per specie' is implemented.

    is_deterministic: bool
        flag to switch between selction criteria

    seed: int
        if is_deterministic==False, seed for the random selection

    """

    def __init__(self, representation, Nselect, act_on='sample per specie', is_deterministic=True, seed=10):
        super(CURFilter, self).__init__()
        self._representation = representation
        self.Nselect = Nselect
        if act_on in ['sample', 'sample per specie', 'feature']:
            self.act_on = act_on
        else:
            raise 'Wrong input: {}'.format(act_on)
        self.is_deterministic = is_deterministic
        self.seed = seed
        self.selected_ids = None

    def fit_transform(self, managers):
        """Perform CUR selection of samples/features.

        Parameters
        ----------
        managers : AtomsList
            list of structures containing features computed with representation

        Returns
        -------
        PseudoPoints
            Selected samples

        Raises
        ------
        ValueError
            [description]
        NotImplementedError
            [description]
        """
        if self.act_on in ['sample per specie']:
            # get the dense feature matrix
            X = managers.get_features(self._representation)

            sps = list(self.Nselect.keys())

            # get various info from the structures about the center atom species and indexing
            strides_by_sp, global_counter, map_by_manager, indices_by_sp =
            self.get_index_mappings_sample_per_species(managers)

            print('The number of pseudo points selected by central atom species is: {}'.format(
                self.Nselect))

            # organize features w.r.t. central atom type
            X_by_sp = {}
            for sp in sps:
                X_by_sp[sp] = X[indices_by_sp[sp]]
            self._XX = X_by_sp

            # split the dense feature matrix by center species and apply CUR decomposition
            selected_ids_by_sp = {}
            for sp in sps:
                print('Selecting species: {}'.format(sp))
                selected_ids_by_sp[sp] = np.sort(do_CUR(X_by_sp[sp], self.Nselect[sp], self.act_on,
                                                        self.is_deterministic, self.seed))

            self.selected_ids = self.convert_selected_global_index2rascal_sample_per_species(
                managers, selected_ids_by_sp, strides_by_sp, map_by_manager)

            # build the pseudo points
            pseudo_points = PseudoPoints(self._representation)
            pseudo_points.extend(managers, self.selected_ids)

            return pseudo_points
        else:
            raise NotImplementedError("method: {}".format(self.act_on))

    def get_index_mappings_sample_per_species(self, managers):
        # get various info from the structures about the center atom species and indexing
        sps = list(self.Nselect.keys())
        types = []
        strides_by_sp = {sp: [0] for sp in sps}
        global_counter = {sp: 0 for sp in sps}
        indices_by_sp = {sp: [] for sp in sps}
        map_by_manager = [{} for ii in range(len(managers))]
        for i_man, man in enumerate(managers):
            counter = {sp: 0 for sp in sps}
            for i_at, at in enumerate(man):
                types.append(at.atom_type)
                if at.atom_type in sps:
                    map_by_manager[i_man][global_counter[at.atom_type]] = i_at
                    counter[at.atom_type] += 1
                    global_counter[at.atom_type] += 1
                else:
                    raise ValueError('Atom type {} has not been specified in fselect: {}'.format(
                        at.atom_type, self.Nselect))
            for sp in sps:
                strides_by_sp[sp].append(counter[sp])

        for sp in sps:
            strides_by_sp[sp] = np.cumsum(strides_by_sp[sp])

        for ii, sp in enumerate(types):
            indices_by_sp[sp].append(ii)

        return strides_by_sp, global_counter, map_by_manager, indices_by_sp

    def convert_selected_global_index2rascal_sample_per_species(self, managers, selected_ids_by_sp, strides_by_sp, map_by_manager):
        # convert selected center indexing into the rascal format
        selected_ids = [[] for ii in range(len(managers))]
        sps = list(self.Nselect.keys())
        i_manager = {sp: 0 for sp in sps}
        for sp in sps:
            for idx in selected_ids_by_sp[sp]:
                carry_on = True
                while carry_on:
                    if idx >= strides_by_sp[sp][i_manager[sp]] and idx < strides_by_sp[sp][i_manager[sp] + 1]:
                        selected_ids[i_manager[sp]].append(
                            map_by_manager[i_manager[sp]][idx])
                        carry_on = False
                    else:
                        i_manager[sp] += 1
        for ii in range(len(selected_ids)):
            selected_ids[ii] = list(np.sort(selected_ids[ii]))
        return selected_ids
