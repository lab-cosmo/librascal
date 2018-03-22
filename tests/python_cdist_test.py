import unittest
import numpy as np
from python_import_proteus import _proteus as pt

class TestCdist(unittest.TestCase):
    def setUp(self):
        """builds the test case. we'll use it to create the matrices
        coordinate matrices A and B between which we wish to compute
        the distances

        """
        self.A = np.array([[0., 0.],
                           [1., 0.],
                           [2., 0.]])
        nb_A = self.A.shape[0]
        self.B = np.array([[0., 1.],
                           [1., 1.]])
        nb_B = self.B.shape[0]

        # the distance matrix is trivial to compute:

        self.dists_ref = np.empty([nb_A, nb_B])

        for i in range(nb_A):
            for j in range(nb_B):
                self.dists_ref[i, j] = np.linalg.norm(
                    self.A[i, :] - self.B[j, :])


    def test_cdist(self):
        """feeds the matrices A and B to Proteus' cdist function and compares
        the results to the local reference dist_ref
        """
        dists = pt.cdist(self.A, self.B)

        error = np.linalg.norm(dists-self.dists_ref)
        tol = 1e-10
        self.assertLessEqual(error, tol)
