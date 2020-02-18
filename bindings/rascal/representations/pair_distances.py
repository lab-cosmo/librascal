import json

from .base import CalculatorFactory, cutoff_function_dict_switch
from ..neighbourlist import AtomsList
import numpy as np


class PairDistances(object):

    """
    Computes a PairDistances representation,

    Hyperparameters
    ----------
    interaction_cutoff : float
        Maximum pairwise distance for atoms to be considered in
        expansion

    cutoff_smooth_width : float
        The distance over which the the interaction is smoothed to zero

    n_species : int
        Number of species to be considered separately

    normalize : boolean
        Whether to normalize so that the kernel between identical environments
        is 1.  Default and highly recommended: True.

    Methods
    -------
    transform(frames)
        Compute the representation for a list of ase.Atoms object.


    .. [soap] Bartók, Kondor, and Csányi, "On representing chemical
        environments", Phys. Rev. B. 87(18), p. 184115
        http://link.aps.org/doi/10.1103/PhysRevB.87.184115

    """

    def __init__(self, interaction_cutoff, cutoff_smooth_width,
                 n_species=1, cutoff_function_type="ShiftedCosine",
                 normalize=True, distance_powers=[],
                 cutoff_function_parameters=dict()):
        """Construct a PairDistances representation

        Required arguments are all the hyperparameters named in the
        class documentation
        """
        self.name = 'pairdistances'
        self.hypers = dict()
        self.update_hyperparameters(
            n_species=n_species, normalize=normalize)
        if distance_powers:
            self.update_hyperparameters(distance_powers)

        cutoff_function_parameters.update(
            interaction_cutoff=interaction_cutoff,
            cutoff_smooth_width=cutoff_smooth_width
        )
        cutoff_function = cutoff_function_dict_switch(cutoff_function_type,
                                **cutoff_function_parameters)

        self.update_hyperparameters(cutoff_function=cutoff_function)

        self.nl_options = [
            dict(name='centers', args=[]),
            dict(name='neighbourlist', args=dict(cutoff=interaction_cutoff)),
            dict(name='strict', args=dict(cutoff=interaction_cutoff))
        ]

        hypers_str = json.dumps(self.hypers)
        self.rep_options = dict(name=self.name, args=[hypers_str])

        self._representation = CalculatorFactory(self.rep_options)

    def update_hyperparameters(self, **hypers):
        """Store the given dict of hyperparameters

        Also updates the internal json-like representation

        """
        allowed_keys = {'interaction_cutoff', 'cutoff_smooth_width',
                        'n_species', 'cutoff_function', 
                        'normalize', 'cutoff_function_parameters'}
        hypers_clean = {key: hypers[key] for key in hypers
                        if key in allowed_keys}
        self.hypers.update(hypers_clean)
        return

    def transform(self, frames):
        """Compute the representation.

        Parameters
        ----------
        frames : list(ase.Atoms)
            List of atomic structures.

        Returns
        -------
            Object containing the representation

        """
        if not isinstance(frames, AtomsList):
            frames = AtomsList(frames, self.nl_options)

        self._representation.compute(frames.managers)

        return frames
