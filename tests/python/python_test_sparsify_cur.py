"""Test the CUR sparsification utilities (plain and wrapped)"""


import ase.io
import collections
import numpy as np
import unittest

from rascal.utils import cur
from rascal import representations


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
        self.assertCountEqual(sel_rows, range(n_sel_rows))
        self.assertCountEqual(sel_cols, range(n_sel_cols))

    def testCURMid(self):
        """Test CUR on a matrix where M < k < N (or vice versa)"""
        n_sel = 450
        sel_rows = cur.do_CUR(self.mat, n_sel, 'sample', verbose=False)
        sel_cols = cur.do_CUR(self.mat.T, n_sel, 'feature', verbose=False)
        self.assertCountEqual(sel_rows, sel_cols)

    def testCURActOn(self):
        """Test behaviour of the 'act_on' argument"""
        n_sel = 200
        sel_rows = cur.do_CUR(self.mat, n_sel, 'rrrsampledsasfd', verbose=False)
        sel_cols = cur.do_CUR(self.mat, n_sel, 'ffffeaturewombat', verbose=False)
        self.assertRaises(ValueError, cur.do_CUR, self.mat, n_sel, 'featuresample')
        self.assertRaises(ValueError, cur.do_CUR, self.mat, n_sel, 'spamspam')

    def testIsDeterministic(self):
        """Test that deterministic selection is independent of the random seed"""
        n_sel = 200
        seed1 = 1618
        sel1 = cur.do_CUR(self.mat, n_sel, 'sample', True, seed1, verbose=False)
        seed2 = 31415
        sel2 = cur.do_CUR(self.mat, n_sel, 'sample', True, seed2, verbose=False)
        self.assertEqual(list(sel1), list(sel2))

    def testFullSparseSVD(self):
        n_sel = 200
        sel1 = cur.do_CUR(self.mat, n_sel, 'sample', True,
                          use_sparse_svd=True, verbose=False)
        sel2 = cur.do_CUR(self.mat, n_sel, 'sample', True,
                          use_sparse_svd=False, verbose=False)
        self.assertCountEqual(sel1, sel2)


class TestCURFilter(unittest.TestCase):
    """Test the CUR wrapper class for rascal representations"""

    def setUp(self):
        self.rep = representations.SphericalInvariants(
            interaction_cutoff=5.0,
            cutoff_smooth_width=0.5,
            max_radial=8,
            max_angular=0,
            gaussian_sigma_type='Constant',
            gaussian_sigma_constant=0.5,
        )
        self.structures = ase.io.read(
            'reference_data/inputs/small_molecules-1000.xyz',
            ':100'
        )
        self.features = self.rep.transform(self.structures)
        self.n_centres = sum(len(struct) for struct in self.structures)

    def testSingleSpecies(self):
        #TODO how do we select the rows of the feature matrix that
        #     correspond to only one type of central atom?
        self.assertTrue(False)

    def testAllSpeciesCentres(self):
        centre_counts = collections.Counter()
        for geom in self.structures:
            centre_counts.update(geom.get_atomic_numbers())
        filter = cur.CURFilter(self.rep, centre_counts,
                               act_on='sample per species')
        filter.filter(self.features)
        n_selected_ids = sum(len(ids) for ids in filter.selected_ids)
        self.assertEqual(n_selected_ids, self.n_centres)

    # TODO test feature selection, zero-of-a-given-species selection,
    #      verify results are the same as returned by plain do_CUR()

