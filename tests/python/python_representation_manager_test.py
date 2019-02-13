import unittest
import numpy as np
import sys
import json

sys.path.insert(0,'../tests/')

from test_utils import load_json_frame, BoxList, Box
import rascal.lib._rascal as rc



class TestSortedCoulombRepresentation(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = '../tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frame = load_json_frame(fn)
        self.structure = unpack_ase(self.frame)
        self.cutoffs = [3.]*self.Natom
        self.cutoff = np.max(self.cutoffs)

        self.nl_options = [
            dict(name='centers',args=[]),
            dict(name='neighbourlist',args=[cutoff]),
            dict(name='strict',args=[cutoff])
        ]

        self.manager =  get_neighbourlist(self.frame,self.nl_options)

        self.central_decay = 0.5
        self.interaction_cutoff = 10.
        self.interaction_decay = 20.
        self.size = 50

        self.inp = json.dumps(
                        dict(central_decay=self.central_decay,
                        interaction_cutoff=self.interaction_cutoff,
                        interaction_decay=self.interaction_decay,
                        size=self.size)
                        )
    def test_manager_iteration(self):
        
