import unittest
import numpy as np
import sys
import json

sys.path.insert(0,'../tests/')

from test_utils import load_json_frame, BoxList, Box
from rascal.representation import SortedCoulombMatrix



class TestSortedCoulombRepresentation(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = '../tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frame = load_json_frame(fn)

        self.hypers = dict(cutoff=3., sorting_algorithm='row_norm',
                        size=50, central_decay=0.5,
                        interaction_cutoff=3, interaction_decay=-1)

    def test_representation_transform(self):

        rep = SortedCoulombMatrix(**self.hypers)

        features = rep.transform([self.frame])

        test = features.get_feature_matrix().T


