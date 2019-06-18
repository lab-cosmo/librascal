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

    def test_hyp2f1(self):
        """
        TEST
        """
        inps = np.random.rand(self.Ntest, 4)

        for ii in range(self.Ntest):
            ref = sps.hyp2f1(*inps[ii, :])
            test = rc.math.hyp2f1(*inps[ii, :])

            self.assertTrue(np.allclose(ref, test, atol=self.atol))
