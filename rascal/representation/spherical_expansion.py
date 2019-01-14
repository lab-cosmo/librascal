import numpy as np
import json

from ..utils import get_strict_neighbourlist
from ..lib import RepresentationManager,FeatureManager
from .base import RepresentationFactory

class SphericalExpansion(object):
    """
    Computes the spherical expansion of the neighbour density [soap]

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


    .. [soap] Bartók, Kondor, and Csányi, "On representing chemical
        environments", Phys. Rev. B. 87(18), p. 184115
        http://link.aps.org/doi/10.1103/PhysRevB.87.184115
    """
    def __init__(self,cutoff,size=10,central_decay=-1,
                    interaction_cutoff=10,interaction_decay=-1):
        self.name = 'sphericalexpansion'
        self.options = 'expansion'
        self.cutoff = cutoff
        self.central_decay = central_decay
        self.interaction_cutoff = interaction_cutoff
        self.interaction_decay = interaction_decay
        self.size = int(size)

    def get_params(self):
        params = dict(name=self.name,sort=self.sort,
                    cutoff=self.cutoff,
                      central_decay=self.central_decay,
                        interaction_cutoff=self.interaction_cutoff,
                        interaction_decay=self.interaction_decay,
                        size=self.size)
        return params

    def transform(self,frames):
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

        managers = list(map(get_strict_neighbourlist,frames,[self.cutoff]*Nframe))

        self.size = self.get_size(managers)

        inp = json.dumps(self.get_params())

        Nfeature = self.get_Nfeature()
        features = FeatureManager.Dense_double(Nfeature,inp)

        cms = map(RepresentationFactory(self.name,self.sort),
                     managers,[inp]*Nframe)

        for cm in cms:
            cm.compute()
            features.append(cm)

        return features

    def get_Nfeature(self):
        return int(self.size*(self.size+1)/2)

    def get_size(self,managers):
        Nneigh = []
        for manager in managers:
            for center in manager:
                Nneigh.append(center.size+1)
        size = int(np.max(Nneigh))
        return size
