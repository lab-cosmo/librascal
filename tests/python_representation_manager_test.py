import unittest
import numpy as np
import sys
import json

sys.path.insert(0,'../tests/')

from test_utils import load_json_frame, BoxList, Box
from python_import_rascal import _rascal as rc


class TestSortedCoulombRepresentation(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        # fn = '../tests/reference_data/simple_cubic_8.json'
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

        self.managerC =  rc.StructureManager.Centers()
        self.managerC.update(np.array(self.positions.T,order='F'),
                       self.numbers.reshape(-1,1),
                       np.array(self.cell.T,order='F'),
                       self.pbc[0].reshape(3,1))

        self.managerLC = rc.Adaptor.NeighbourList_Centers(self.managerC,self.max_cutoff)
        self.managerLC.update()

        self.manager = rc.Adaptor.Strict_NeighbourList_Centers(self.managerLC,self.max_cutoff)
        self.manager.update()

        self.central_decay = 0.5
        self.interaction_cutoff = 10.
        self.interaction_decay = 20.
        self.size = 50

        self.inp = json.dumps(
                        dict(central_decay=self.central_decay,
                        interaction_cutoff=self.interaction_cutoff,
                        interaction_decay=self.interaction_decay,
                        size=self.size)
                        )

    def test_constructor1(self):
        """
        TEST constructor wrapper
        """
        rc.RepresentationManager.SortedCoulomb_Strict_NeighbourList_Centers_SortDistance(
                        self.manager,self.inp)

    def test_compute1(self):
        """
        TEST compute wrapper
        """
        cm = rc.RepresentationManager.SortedCoulomb_Strict_NeighbourList_Centers_SortDistance(
                        self.manager,self.inp)
        cm.compute()

    def test_constructor2(self):
        """
        TEST constructor wrapper
        """
        rc.RepresentationManager.SortedCoulomb_Strict_NeighbourList_Centers_SortRowNorm(
                        self.manager,self.inp)

    def test_compute2(self):
        """
        TEST compute wrapper
        """
        cm = rc.RepresentationManager.SortedCoulomb_Strict_NeighbourList_Centers_SortRowNorm(
                        self.manager,self.inp)
        cm.compute()