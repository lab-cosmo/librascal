import numpy as np
import json

from ..neighbourlist import get_neighbourlist, get_neighbourlist_full_name
from ..lib import RepresentationManager,FeatureManager
from .base import RepresentationFactory


class SortedCoulombMatrix(object):
    """
    Computes the Sorted Coulomb matrix representation [1].

    Attributes
    ----------
    cutoff : float

    central_decay : float
        The distance over which the the coulomb interaction decays from full to none.

    interaction_cutoff : float
        The distance between two non-central atom, where the coulomb interaction element will be zero.

    interaction_decay : float
        The distance over which the the coulomb interaction decays from full to none.

    size : int
        Larger or equal to the maximum number of neighbour an atom has in the structure.

    Methods
    -------
    transform(frames)
        Compute the representation for a list of ase.Atoms object.


    .. [1] Rupp, M., Tkatchenko, A., Müller, K.-R., & von Lilienfeld, O. A. (2011).
        Fast and Accurate Modeling of Molecular Atomization Energies with Machine Learning.
        Physical Review Letters, 108(5), 58301. https://doi.org/10.1103/PhysRevLett.108.058301
    """
    def __init__(self, cutoff, sorting_algorithm='row_norm', size=10, central_decay=-1, interaction_cutoff=10, interaction_decay=-1):
        self.name = 'sortedcoulomb'
        self.sorting_algorithm = sorting_algorithm
        self.cutoff = cutoff
        self.central_decay = central_decay
        self.interaction_cutoff = interaction_cutoff
        self.interaction_decay = interaction_decay
        self.size = int(size)

        self.nl_options = [
                dict(name='centers',args=[]),
                dict(name='neighbourlist',args=[cutoff]),
                dict(name='strict',args=[cutoff])
        ]

        neighbourlist_full_name = get_neighbourlist_full_name(self.nl_options)
        self.name = self.name + '_' + neighbourlist_full_name

    def get_params(self):
        params = dict(name=self.name,sorting_algorithm=self.sorting_algorithm,
                    cutoff=self.cutoff,nl_options=self.nl_options,
                    central_decay=self.central_decay,
                    interaction_cutoff=self.interaction_cutoff,
                    interaction_decay=self.interaction_decay,
                    size=self.size)
        return params

    def transform(self, frames):
        """
        Compute the representation.

        Parameters
        ----------
        frames : list(ase.Atoms)
            List of atomic structures.

        Returns
        -------
        FeatureManager.Dense_double
            Object containing the representation
        """
        Nframe = len(frames)

        managers = list(map(get_neighbourlist,frames,[self.nl_options]*Nframe))

        self.size = self.get_size(managers)

        inp = json.dumps(self.get_params())

        Nfeature = self.get_Nfeature()
        features = FeatureManager.Dense_double(Nfeature, inp)

        cms = map(RepresentationFactory,[self.name]*Nframe,
                     managers,[inp]*Nframe)

        for cm in cms:
            cm.compute()
            features.append(cm)

        return features

    def get_Nfeature(self):
        return int(self.size*(self.size+1)/2)

    def get_size(self, managers):
        Nneigh = []
        for manager in managers:
            for center in manager:
                Nneigh.append(center.size+1)
        size = int(np.max(Nneigh))
        return size
