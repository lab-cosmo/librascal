import numpy as np
import unittest

from rascal.utils import cur


class TestCURFun(unittest.TestCase):

    def setUp(self):
        np.random.seed(2718)
        self.mat = np.random.random((400, 500))

    def testCURNormal(self):
        """Test CUR on a matrix: Usual case (k < N and k < M)"""
        n_sel = 200
        sel_rows = cur.do_CUR(self.mat, n_sel, 'sample')
        sel_cols = cur.do_CUR(self.mat, n_sel, 'feature')
        self.assertEqual(len(sel_rows), n_sel, "Wrong length of selection")
        self.assertEqual(len(sel_cols), n_sel, "Wrong length of selection")

    def testCURLarge(self):
        """Test CUR on a matrix where k = M or N"""
        n_sel_rows = 400
        n_sel_cols = 500
        sel_rows = cur.do_CUR(self.mat, n_sel_rows, 'sample')
        sel_cols = cur.do_CUR(self.mat, n_sel_cols, 'feature')
        self.assertEqual(sel_rows, np.arange(n_sel_rows))
        self.assertEqual(sel_cols, np.arange(n_sel_cols))

    def testCURMid(self):
        """Test CUR on a matrix where M < k < N (or vice versa)"""
        n_sel = 450
        sel_cols = cur.do_CUR(self.mat, n_sel, 'feature')
        sel_rows = cur.do_CUR(self.mat.T, n_sel, 'sample')
        self.assertEqual(sel_rows, sel_cols)

    def testCURActOn(self):
        """Test behaviour of the 'act_on' argument"""
        n_sel = 200
        sel_rows = cur.do_CUR(self.mat, n_sel, 'rrrsampledsasfd')
        sel_cols = cur.do_CUR(self.mat, n_sel, 'ffffeaturewombat')
        self.assertRaises(ValueError, cur.do_CUR, self.mat, n_sel, 'featuresample')
        self.assertRaises(ValueError, cur.do_CUR, self.mat, n_sel, 'blablafoobar')
        return

    def testIsDeterministic(self):
        """Test that deterministic selection is independent of the random seed"""
        n_sel = 200
        seed1 = 1618
        sel1 = cur.do_CUR(self.mat, n_sel, 'sample', True, seed1)
        seed2 = 1618
        sel2 = cur.do_CUR(self.mat, n_sel, 'sample', True, seed2)
        self.assertEqual(sel1, sel2)
        return

    # TODO test full vs. sparse SVD (once implemented)

