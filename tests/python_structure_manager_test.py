import unittest
import numpy as np
import sys
import faulthandler

sys.path.insert(0,'../tests/')

from test_utils import load_json_frame, BoxList, Box
from python_import_rascal import _rascal as rc


def get_NL_reference(cutoff,cell,pbc,positions,numbers):
    list_box = BoxList(cutoff,cell,pbc,positions)

    neighlist = [[] for it in range(len(numbers))]
    neighpos = [[] for it in range(len(numbers))]
    neighshift = [[] for it in range(len(numbers))]
    neighdist = [[] for it in range(len(numbers))]
    neightype = [[] for it in range(len(numbers))]
    dirVec = [[] for it in range(len(numbers))]
    for box in list_box.iter_box():
        for icenter in box.icenters:
            for jneigh,box_shift in box.iter_neigh_box():
            
                nnp = positions[jneigh]+ np.dot(box_shift.reshape((1,3)),cell).reshape((1,3))
                rr = nnp - positions[icenter].reshape((1,3))
                dist = np.linalg.norm(rr)

                neighpos[icenter].extend(nnp)
                neighlist[icenter].append(jneigh)
                neightype[icenter].append(numbers[jneigh])
                neighdist[icenter].append(dist)
                
    return neighpos,neighlist,neightype,neighdist

def get_NL_strict_reference(cutoff,cell,pbc,positions,numbers):
    list_box = BoxList(cutoff,cell,pbc,positions)

    neighlist = [[] for it in range(len(numbers))]
    neighpos = [[] for it in range(len(numbers))]
    neighshift = [[] for it in range(len(numbers))]
    neighdist = [[] for it in range(len(numbers))]
    neightype = [[] for it in range(len(numbers))]
    dirVec = [[] for it in range(len(numbers))]
    for box in list_box.iter_box():
        for icenter in box.icenters:
            for jneigh,box_shift in box.iter_neigh_box():
            
                nnp = positions[jneigh]+ np.dot(box_shift.reshape((1,3)),cell).reshape((1,3))
                rr = nnp - positions[icenter].reshape((1,3))
                dist = np.linalg.norm(rr)

                if cutoff > dist and dist > 0:
                    neighpos[icenter].extend(nnp)
                    neighlist[icenter].append(jneigh)
                    neighdist[icenter].append(dist)
                    neightype[icenter].append(numbers[jneigh])
                    dirVec[icenter].append(rr/dist)
    return neighpos,neighlist,neightype,neighdist,dirVec


class TestStructureManagerCenters(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = '../tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frame = load_json_frame(fn)

        self.cell = self.frame['cell']
        self.positions = self.frame['positions']
        self.numbers = self.frame['numbers']
        self.pbc = np.array([ [1, 1, 1], [0, 0, 0],
                              [0, 1, 0], [1, 0, 1],
                              [1, 1, 0], [0, 0, 1],
                              [1, 0, 0], [0, 1, 0] ]).astype(int)

        self.Natom = self.positions.shape[0]
        self.cutoffs = [3.]*self.Natom
        self.max_cutoff = np.max(self.cutoffs)

    def test_constructor(self):
        """
        TEST constructor wrapper
        """
        rc.StructureManager.Centers()
    
    def test_update(self):
        """
        TEST constructor wrapper
        """
        manager =  rc.StructureManager.Centers()
        centers = np.array([it for it in range(self.Natom)], dtype=np.int32)
        manager.update(np.array(self.positions.T,order='F'),
                       self.numbers.reshape(-1,1),
                       np.array(self.cell.T,order='F'),
                       self.pbc[0].reshape(3,1))

    def test_manager_iteration(self):
        manager =  rc.StructureManager.Centers()
        centers = np.array([it for it in range(self.Natom)], dtype=np.int32)
        manager.update(np.array(self.positions.T,order='F'),
                       self.numbers.reshape(-1,1),
                       np.array(self.cell.T,order='F'),
                       self.pbc[0].reshape(3,1))
        ii = 0
        for center in manager:
            self.assertTrue(ii == center.atom_index)
            self.assertTrue(self.numbers[ii] == center.atom_type)
            self.assertTrue(np.allclose(self.positions[ii], center.position))
            ii += 1

class TestNL(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = '../tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frame = load_json_frame(fn)

        self.cell = self.frame['cell']
        self.positions = self.frame['positions']
        self.numbers = self.frame['numbers']
        self.pbcs = np.array([ [1, 1, 1], [0, 0, 0],
                              [0, 1, 0], [1, 0, 1],
                              [1, 1, 0], [0, 0, 1],
                              [1, 0, 0], [0, 1, 0] ]).astype(int)

        self.Natom = self.positions.shape[0]
        self.cutoffs = [3.]*self.Natom
        self.max_cutoff = np.max(self.cutoffs)

        self.manager =  rc.StructureManager.Centers()
        self.manager.update(np.array(self.positions.T,order='F'),
                       self.numbers.reshape(-1,1),
                       np.array(self.cell.T,order='F'),
                       self.pbcs[0].reshape(3,1))

    def test_constructor(self):
        """
        TEST constructor wrapper
        """
        rc.Adaptor.NeighbourList_Centers(self.manager,self.max_cutoff)
    
    def test_update(self):
        manager =  rc.Adaptor.NeighbourList_Centers(self.manager,self.max_cutoff)
        manager.update()

    def test_manager_iteration_1(self):
        manager =  rc.Adaptor.NeighbourList_Centers(self.manager,self.max_cutoff)
        manager.update()

        ii = 0
        for center in manager:
            self.assertTrue(ii == center.atom_index)
            self.assertTrue(self.numbers[ii] == center.atom_type)
            self.assertTrue(np.allclose(self.positions[ii], center.position))
            ii += 1

    def test_manager_iteration_2(self):
        manager =  rc.Adaptor.NeighbourList_Centers(self.manager,self.max_cutoff)
        manager.update()
        
        neighpos,neighlist,neightype,neighdist = get_NL_reference(
                    self.max_cutoff,self.cell,self.pbcs[0],self.positions,self.numbers)

        for ii,center in enumerate(manager):
            for jj,neigh in enumerate(center):
                dist = np.linalg.norm(neigh.position-center.position)

class TestNLStrict(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = '../tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frame = load_json_frame(fn)

        self.cell = self.frame['cell']
        self.positions = self.frame['positions']
        self.numbers = self.frame['numbers']
        self.pbcs = np.array([ [1, 1, 1], [0, 0, 0],
                              [0, 1, 0], [1, 0, 1],
                              [1, 1, 0], [0, 0, 1],
                              [1, 0, 0], [0, 1, 0] ]).astype(int)

        self.Natom = self.positions.shape[0]
        self.cutoffs = [3.]*self.Natom
        self.max_cutoff = np.max(self.cutoffs)

        self.managerC =  rc.StructureManager.Centers()
        self.managerC.update(np.array(self.positions.T,order='F'),
                       self.numbers.reshape(-1,1),
                       np.array(self.cell.T,order='F'),
                       self.pbcs[0].reshape(3,1))
        self.manager = rc.Adaptor.NeighbourList_Centers(self.managerC,self.max_cutoff) 
        self.manager.update()

    def test_constructor(self):
        """
        TEST constructor wrapper
        """
        rc.Adaptor.Strict_NeighbourList_Centers(self.manager,self.max_cutoff)
    
    def test_a_update(self):
        manager =  rc.Adaptor.Strict_NeighbourList_Centers(self.manager,self.max_cutoff)
        manager.update()

    def test_manager_iteration_1(self):
        manager = rc.Adaptor.Strict_NeighbourList_Centers(self.manager,self.max_cutoff)
        manager.update()

        ii = 0
        for center in manager:
            self.assertEqual(ii,center.atom_index)
            self.assertTrue(self.numbers[ii] == center.atom_type)
            self.assertTrue(np.allclose(self.positions[ii], center.position))
            ii += 1
    
    def test_manager_iteration_2(self):
        '''
        Compare the distances and direction vector between the reference and 
        librascal sctrict neighbourlist 
        '''
        manager = rc.Adaptor.Strict_NeighbourList_Centers(self.manager,self.max_cutoff)
        manager.update()

        neighpos,neighlist,neightype,neighdist,neighdirVec = get_NL_strict_reference(
                    self.max_cutoff,self.cell,self.pbcs[0],self.positions,self.numbers)

        for ii,center in enumerate(manager):
            dists,dirVecs = [],[]
            for jj,neigh in enumerate(center):
                dist = np.linalg.norm(neigh.position-center.position)
                dists.append(dist)
                dirVecs.append((neigh.position-center.position)/dist)

            ref_dists, dists = np.array(neighdist[ii]),np.array(dists)
            ref_dirVecs, dirVecs = np.array(neighdirVec[ii]).reshape((-1,3)),np.array(dirVecs)
            # sort because the order is not the same
            ref_sort_ids,sort_ids = np.argsort(ref_dists),np.argsort(dists)
            self.assertTrue(np.allclose(ref_dists[ref_sort_ids],dists[sort_ids]))
            self.assertTrue(np.allclose(ref_dirVecs[ref_sort_ids],dirVecs[sort_ids]))
            
