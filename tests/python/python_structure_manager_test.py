from rascal.neighbourlist import get_neighbourlist
from test_utils import load_json_frame, BoxList, Box
import unittest
import numpy as np
import sys
import faulthandler


def get_NL_reference(cutoff, cell, pbc, positions, atom_types):
    list_box = BoxList(cutoff, np.array(cell.T, order='C'),
                       pbc.flatten(), np.array(positions.T, order='C'))

    neighlist = [[] for it in range(len(atom_types))]
    neighpos = [[] for it in range(len(atom_types))]
    neighshift = [[] for it in range(len(atom_types))]
    neighdist = [[] for it in range(len(atom_types))]
    neightype = [[] for it in range(len(atom_types))]
    dirVec = [[] for it in range(len(atom_types))]
    for box in list_box.iter_box():
        for icenter in box.icenters:
            for jneigh, box_shift in box.iter_neigh_box():

                nnp = positions[:, jneigh] + \
                    np.dot(box_shift.reshape((1, 3)), cell).reshape((1, 3))
                rr = nnp - positions[:, icenter].reshape((1, 3))
                dist = np.linalg.norm(rr)

                neighpos[icenter].extend(nnp)
                neighlist[icenter].append(jneigh)
                neightype[icenter].append(atom_types[jneigh])
                neighdist[icenter].append(dist)

    return neighpos, neighlist, neightype, neighdist


def get_NL_strict_reference(cutoff, cell, pbc, positions, atom_types):
    list_box = BoxList(cutoff, np.array(cell.T, order='C'),
                       pbc.flatten(), np.array(positions.T, order='C'))

    neighlist = [[] for it in range(len(atom_types))]
    neighpos = [[] for it in range(len(atom_types))]
    neighshift = [[] for it in range(len(atom_types))]
    neighdist = [[] for it in range(len(atom_types))]
    neightype = [[] for it in range(len(atom_types))]
    dirVec = [[] for it in range(len(atom_types))]
    for box in list_box.iter_box():
        for icenter in box.icenters:
            for jneigh, box_shift in box.iter_neigh_box():

                nnp = positions[:, jneigh] + \
                    np.dot(box_shift.reshape((1, 3)), cell).reshape((1, 3))
                rr = nnp - positions[:, icenter].reshape((1, 3))
                dist = np.linalg.norm(rr)

                if cutoff > dist and dist > 0:
                    neighpos[icenter].extend(nnp)
                    neighlist[icenter].append(jneigh)
                    neighdist[icenter].append(dist)
                    neightype[icenter].append(atom_types[jneigh])
                    dirVec[icenter].append(rr/dist)
    return neighpos, neighlist, neightype, neighdist, dirVec


class TestStructureManagerCenters(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = '../tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frame = load_json_frame(fn)
        self.structure = self.frame
        self.nl_options = [
            dict(name='centers', args={}),
        ]

    def test_manager_iteration(self):
        manager = get_neighbourlist(self.frame, self.nl_options)
        ii = 0
        for center in manager:
            self.assertTrue(ii == center.atom_tag)
            self.assertTrue(
                self.structure['atom_types'][ii] == center.atom_type)
            self.assertTrue(np.allclose(
                self.structure['positions'][:, ii], center.position))
            ii += 1


class TestNL(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = '../tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frame = load_json_frame(fn)
        self.structure = self.frame
        self.cutoff = 3.
        self.nl_options = [
            dict(name='centers', args=dict()),
            dict(name='neighbourlist', args=dict(cutoff=self.cutoff)),
        ]
        self.pbcs = np.array([[1, 1, 1], [0, 0, 0],
                              [0, 1, 0], [1, 0, 1],
                              [1, 1, 0], [0, 0, 1],
                              [1, 0, 0], [0, 1, 0]]).astype(int)

    def test_manager_iteration_1(self):
        manager = get_neighbourlist(self.frame, self.nl_options)
        ii = 0
        for center in manager:
            self.assertTrue(ii == center.atom_tag)
            self.assertTrue(
                self.structure['atom_types'][ii] == center.atom_type)
            self.assertTrue(np.allclose(
                self.structure['positions'][:, ii], center.position))
            ii += 1

    def test_manager_iteration_2(self):
        frame = self.frame
        structure = self.structure
        for pbc in self.pbcs:
            frame['pbc'] = pbc
            structure['pbc'] = pbc
            manager = get_neighbourlist(frame, self.nl_options)

            neighpos, neighlist, neightype, neighdist = get_NL_reference(
                self.cutoff, **structure)

            for ii, center in enumerate(manager):
                for jj, neigh in enumerate(center):
                    dist = np.linalg.norm(neigh.position-center.position)


class TestNLStrict(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = '../tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frame = load_json_frame(fn)
        self.structure = self.frame
        self.cutoff = 3.

        self.nl_options = [
            dict(name='centers', args=dict()),
            dict(name='neighbourlist', args=dict(cutoff=self.cutoff)),
            dict(name='strict', args=dict(cutoff=self.cutoff))
        ]

        self.pbcs = np.array([[1, 1, 1], [0, 0, 0],
                              [0, 1, 0], [1, 0, 1],
                              [1, 1, 0], [0, 0, 1],
                              [1, 0, 0], [0, 1, 0]]).astype(int)

    def test_manager_iteration_1(self):
        manager = get_neighbourlist(self.frame, self.nl_options)
        ii = 0
        for center in manager:
            self.assertTrue(ii == center.atom_tag)
            self.assertTrue(
                self.structure['atom_types'][ii] == center.atom_type)
            self.assertTrue(np.allclose(
                self.structure['positions'][:, ii], center.position))
            ii += 1

    def test_manager_iteration_2(self):
        '''
        Compare the distances and direction vector between the reference and
        librascal sctrict neighbourlist
        '''
        frame = self.frame
        structure = self.structure
        for pbc in self.pbcs:
            frame['pbc'] = pbc
            structure['pbc'] = pbc
            manager = get_neighbourlist(frame, self.nl_options)

            neighpos, neighlist, neightype, neighdist, neighdirVec = get_NL_strict_reference(
                self.cutoff, **structure)
            for ii, center in enumerate(manager):
                dists, dirVecs = [], []
                for jj, neigh in enumerate(center):
                    dist = np.linalg.norm(neigh.position-center.position)
                    dists.append(dist)
                    dirVecs.append((neigh.position-center.position)/dist)

                ref_dists, dists = np.array(neighdist[ii]), np.array(dists)
                ref_dirVecs, dirVecs = np.array(
                    neighdirVec[ii]).reshape((-1, 3)), np.array(dirVecs)
                # sort because the order is not the same
                ref_sort_ids, sort_ids = np.argsort(
                    ref_dists), np.argsort(dists)
