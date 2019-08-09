import rascal.lib._rascal as rc
import unittest
import numpy as np
import scipy.special as sps
import sys


class TestMath(unittest.TestCase):
    def setUp(self):
        """

        """
        np.random.seed(10)
        self.Ntest = 100
        self.atol = 1e-14
