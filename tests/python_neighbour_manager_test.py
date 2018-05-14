import unittest
import numpy as np
import sys
sys.path.insert(0,'../tests/')
from test_utils import load_json_frame,BoxList,Box
from python_import_rascal import _rascal as  rc

class TestNeighbourManagerCell(unittest.TestCase):
    def setUp(self):
        """builds the test case. Test the cell neighbourlist implementation against a triclinic crystal using 
        a custom pyhton implementation as a reference (which was tested against ase neighbourlist).

        """

        fn = '../tests/reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frame = load_json_frame(fn)
        
        self.cell = self.frame['cell']
        self.positions = self.frame['positions']
        self.numbers = self.frame['numbers']
        self.pbc = [[True,True,True],[False,False,False],
                    [False,True,True],[True,False,True],[True,True,False],
                    [False,False,True],[True,False,False],[False,True,False]]
        self.Natom = self.positions.shape[0]
        self.cutoffs = [3.]*self.Natom
        self.max_cutoff = np.max(self.cutoffs)

    def get_reference(self,pbc):
        list_box = BoxList(self.max_cutoff,self.cell,pbc,self.positions)
        neighlist = [[] for it in range(self.Natom)]
        neighpos = [[] for it in range(self.Natom)]
        neighshift = [[] for it in range(self.Natom)]
        for box in list_box.iter_box():
            for icenter in box.icenters:
                for jneigh,box_shift in box.iter_neigh_box():
                
                    nnp = self.positions[jneigh] + np.dot(box_shift.reshape((1,3)),self.cell).reshape((1,3))
                    #rr = nnp - self.positions[icenter].reshape((1,3))
                    #dist = np.linalg.norm(rr,axis=1)
                    neighpos[icenter].extend(nnp)
                    neighlist[icenter].append(jneigh)
                    neighshift[icenter].append(box_shift)
        return neighpos,neighlist,neighshift

    def test_constructor(self):
        """
        TEST constructor wrapper
        """
        rc.NeighbourhoodManagerCell()
    
    def test_manager_iteration(self):
        manager =  rc.NeighbourhoodManagerCell()
        centers = np.array([it for it in range(self.Natom)],dtype=np.int32)
        manager.update(np.array(self.positions.T,order='F'),self.numbers,centers, 
                      np.array(self.cell.T,order='F'),self.pbc[0],self.max_cutoff)
        ii = 0
        for center in manager:
            self.assertTrue(ii == center.get_atom_index())
            self.assertTrue(self.numbers[ii] == center.get_atom_type())
            self.assertTrue(np.allclose(self.positions[ii], center.get_position()))
            ii += 1
        
    def test_neighbour_iteration(self):
        for pbc in self.pbc:
            neighpos,neighlist,neighshift = self.get_reference(pbc)

            manager =  rc.NeighbourhoodManagerCell()
            centers = np.array([it for it in range(self.Natom)],dtype=np.int32)
            manager.update(np.array(self.positions.T,order='F'),self.numbers,centers, 
                            np.array(self.cell.T,order='F'),pbc,self.max_cutoff)
            
            for center in manager:
                icenter = center.get_atom_index()
                ii = 0
                for neigh in center:
                    ineigh = neigh.get_atom_index()
                    self.assertTrue(neighlist[icenter][ii] == neigh.get_atom_index())
                    self.assertTrue(self.numbers[ineigh] == neigh.get_atom_type())
                    
                    self.assertTrue(np.allclose(neighpos[icenter][ii], neigh.get_position()))
                    ii += 1




