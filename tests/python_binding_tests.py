#!/usr/bin/env python3

import unittest

from python_cdist_test import TestCdist
from python_neighbour_cell_test import TestNeighbourManagerCell

class SimpleCheck(unittest.TestCase):
    def setUp(self):
        self.truth = True


    def test_simple_example(self):
        self.assertTrue(self.truth)


if __name__ == '__main__':
    unittest.main()

