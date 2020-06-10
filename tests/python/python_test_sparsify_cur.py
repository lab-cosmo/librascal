"""Test the CUR sparsification utilities (plain and wrapped)"""


import numpy as np
import unittest

from rascal.utils import cur


class TestCURFun(unittest.TestCase):
    """Test the CUR function for plain matrices ('do_CUR')"""

    def setUp(self):
        np.random.seed(2718)
        self.mat = np.random.random((400, 500))

    def testCURNormal(self):
        """Test CUR on a matrix: Usual case (k < N and k < M)"""
        n_sel = 200
        sel_rows = cur.do_CUR(self.mat, n_sel, 'sample', verbose=False)
        sel_cols = cur.do_CUR(self.mat, n_sel, 'feature', verbose=False)
        self.assertEqual(len(sel_rows), n_sel, "Wrong length of selection")
        self.assertEqual(len(sel_cols), n_sel, "Wrong length of selection")

    def testCURLarge(self):
        """Test CUR on a matrix where k = M or N"""
        n_sel_rows = 400
        n_sel_cols = 500
        sel_rows = cur.do_CUR(self.mat, n_sel_rows, 'sample', verbose=False)
        sel_cols = cur.do_CUR(self.mat, n_sel_cols, 'feature', verbose=False)
        self.assertEqual(sorted(list(sel_rows)), list(range(n_sel_rows)))
        self.assertEqual(sorted(list(sel_cols)), list(range(n_sel_cols)))

    def testCURMid(self):
        """Test CUR on a matrix where M < k < N (or vice versa)"""
        n_sel = 450
        sel_rows = cur.do_CUR(self.mat, n_sel, 'sample', verbose=False)
        sel_cols = cur.do_CUR(self.mat.T, n_sel, 'feature', verbose=False)
        self.assertEqual(sorted(list(sel_rows)), sorted(list(sel_cols)))

    def testCURActOn(self):
        """Test behaviour of the 'act_on' argument"""
        n_sel = 200
        sel_rows = cur.do_CUR(self.mat, n_sel, 'rrrsampledsasfd', verbose=False)
        sel_cols = cur.do_CUR(self.mat, n_sel, 'ffffeaturewombat', verbose=False)
        self.assertRaises(ValueError, cur.do_CUR, self.mat, n_sel, 'featuresample')
        self.assertRaises(ValueError, cur.do_CUR, self.mat, n_sel, 'blablafoobar')

    def testIsDeterministic(self):
        """Test that deterministic selection is independent of the random seed"""
        n_sel = 200
        seed1 = 1618
        sel1 = cur.do_CUR(self.mat, n_sel, 'sample', True, seed1, verbose=False)
        seed2 = 31415
        sel2 = cur.do_CUR(self.mat, n_sel, 'sample', True, seed2, verbose=False)
        self.assertEqual(list(sel1), list(sel2))

    # TODO test full vs. sparse SVD (once implemented)

