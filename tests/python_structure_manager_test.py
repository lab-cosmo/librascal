import unittest
import numpy as np
import sys
import faulthandler

sys.path.insert(0,'../tests/')

from test_utils import load_json_frame, BoxList, Box
from python_import_rascal import _rascal as rc


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
        self.pbc = np.array([ [1, 1, 1], [0, 0, 0],
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
                       self.pbc[0].reshape(3,1))

    def test_constructor(self):
        """
        TEST constructor wrapper
        """
        rc.Adaptor.NeighbourList_Centers(self.manager,self.max_cutoff)
    
    def test_update(self):
        manager =  rc.Adaptor.NeighbourList_Centers(self.manager,self.max_cutoff)
        manager.update()

    def test_manager_iteration(self):
        manager =  rc.Adaptor.NeighbourList_Centers(self.manager,self.max_cutoff)
        manager.update()

        ii = 0
        for center in manager:
            self.assertTrue(ii == center.atom_index)
            self.assertTrue(self.numbers[ii] == center.atom_type)
            self.assertTrue(np.allclose(self.positions[ii], center.position))
            ii += 1

class TestNLStrict(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = '../tests/reference_data/simple_cubic_8.json'
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

        self.managerC =  rc.StructureManager.Centers()
        self.managerC.update(np.array(self.positions.T,order='F'),
                       self.numbers.reshape(-1,1),
                       np.array(self.cell.T,order='F'),
                       self.pbc[0].reshape(3,1))
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
            self.assertTrue(ii == center.atom_index)
            self.assertTrue(self.numbers[ii] == center.atom_type)
            self.assertTrue(np.allclose(self.positions[ii], center.position))
            ii += 1
    
    def test_manager_iteration_2(self):
        manager = rc.Adaptor.Strict_NeighbourList_Centers(self.manager,self.max_cutoff)
        manager.update()

        ii = 0
        for center in manager:
            for neigh in center:
                pass



# class TestStructureManagerCell(unittest.TestCase):
#     def setUp(self):
#         """builds the test case. Test the cell neighbourlist implementation against a triclinic crystal using
#         a custom pyhton implementation as a reference (which was tested against ase neighbourlist).

#         """

#         fn = '../tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
#         self.frame = load_json_frame(fn)

#         self.cell = self.frame['cell']
#         self.positions = self.frame['positions']
#         self.numbers = self.frame['numbers']
#         self.pbc = [[True,True,True],[False,False,False],
#                     [False,True,True],[True,False,True],[True,True,False],
#                     [False,False,True],[True,False,False],[False,True,False]]
#         self.Natom = self.positions.shape[0]
#         self.cutoffs = [3.]*self.Natom
#         self.max_cutoff = np.max(self.cutoffs)

#     def get_reference(self,pbc):
#         list_box = BoxList(self.max_cutoff,self.cell,pbc,self.positions)
#         neighlist = [[] for it in range(self.Natom)]
#         neighpos = [[] for it in range(self.Natom)]
#         neighshift = [[] for it in range(self.Natom)]
#         for box in list_box.iter_box():
#             for icenter in box.icenters:
#                 for jneigh,box_shift in box.iter_neigh_box():

#                     nnp = self.positions[jneigh] + np.dot(box_shift.reshape((1,3)),self.cell).reshape((1,3))
#                     #rr = nnp - self.positions[icenter].reshape((1,3))
#                     #dist = np.linalg.norm(rr,axis=1)
#                     neighpos[icenter].extend(nnp)
#                     neighlist[icenter].append(jneigh)
#                     neighshift[icenter].append(box_shift)
#         return neighpos,neighlist,neighshift

#     def test_constructor(self):
#         """
#         TEST constructor wrapper
#         """
#         rc.StructureManagerCell()

#     def test_manager_iteration(self):
#         manager =  rc.StructureManagerCell()
#         centers = np.array([it for it in range(self.Natom)],dtype=np.int32)
#         manager.update(np.array(self.positions.T,order='F'),self.numbers,centers,
#                       np.array(self.cell.T,order='F'),self.pbc[0],self.max_cutoff)
#         ii = 0
#         for center in manager:
#             self.assertTrue(ii == center.atom_index)
#             self.assertTrue(self.numbers[ii] == center.atom_type)
#             self.assertTrue(np.allclose(self.positions[ii], center.position))
#             ii += 1

#     def test_neighbour_iteration(self):
#         for pbc in self.pbc:
#             neighpos,neighlist,neighshift = self.get_reference(pbc)

#             manager =  rc.StructureManagerCell()
#             centers = np.array([it for it in range(self.Natom)],dtype=np.int32)
#             manager.update(np.array(self.positions.T,order='F'),self.numbers,centers,
#                             np.array(self.cell.T,order='F'),pbc,self.max_cutoff)

#             for center in manager:
#                 icenter = center.atom_index
#                 for ii,neigh in enumerate(center):
#                     ineigh = neigh.atom_index
#                     try:
#                         self.assertTrue(neighlist[icenter][ii] == neigh.atom_index)
#                     except Exception as error:
#                         raise Exception ("neigh.atom_index ={} {}".format(neigh.atom_index, error))



#                     self.assertEqual(self.numbers[ineigh], neigh.atom_type)

#                     val = (neighpos[icenter][ii] == neigh.position).all()
#                     ## TODO: this test fails due to a referenced bug
#                     ## in get_neighbour_position in
#                     ## structure_manager_cell.hh:128
#                     if not val:
#                         #print ("neighpos {}, neigh.pos{}".format(neighpos[icenter][ii], neigh.position))
#                         pass

#                     #self.assertTrue(val)
