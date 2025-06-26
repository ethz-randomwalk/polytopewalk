import unittest
from polytopewalk import FacialReduction
import numpy as np 
import pandas as pd
from scipy.sparse import lil_matrix
from polytopewalk.sparse import *
from polytopewalk.dense import *

simplex_A = lil_matrix((1, 3))
simplex_A[(0, 0)] = 1
simplex_A[(0, 1)] = 1
simplex_A[(0, 2)] = 1
simplex_b = np.array([1])

hc_A = lil_matrix((4,6))
hc_A[(0, 0)] = 1
hc_A[(0, 2)] = 1
hc_A[(1, 1)] = 1
hc_A[(1, 3)] = 1
hc_A[(2, 0)] = -1
hc_A[(2, 4)] = 1
hc_A[(3, 1)] = -1
hc_A[(3, 5)] = 1

hc_b = np.array([1,1,1,1])

birk_A = lil_matrix((3, 4))
birk_A[(0, 0)] = 1
birk_A[(0, 1)] = 1
birk_A[(1, 2)] = 1
birk_A[(1, 3)] = 1
birk_A[(2, 0)] = 1
birk_A[(2, 2)] = 1

birk_b = np.array([1,1,1])

class TestWeights(unittest.TestCase):
    def test_weights(self):
        sparse_vaidya = SparseVaidyaWalk(r = 0.9)
        sparse_john = SparseJohnWalk(r = 0.9)
        sparse_dikinls = SparseDikinLSWalk(r = 0.9)

        simplex_start = np.array([0.33, 0.34, 0.33])
        w = sparse_dikinls.generateWeight(simplex_start, simplex_A, 3)
        self.assertAlmostEqual(w.sum(), 2, places = 1)
        w = sparse_john.generateWeight(simplex_start, simplex_A, 3)
        self.assertAlmostEqual(w.sum(), 3, places = 2)
        w = sparse_vaidya.generateWeight(simplex_start, simplex_A, 3)
        self.assertAlmostEqual(w.sum(), 4, places = 2)

        hc_start = np.array([0, 0, 1, 1, 1, 1])
        w = sparse_dikinls.generateWeight(hc_start, hc_A, 4)
        self.assertAlmostEqual(w.sum(), 2, places = 1)
        w = sparse_john.generateWeight(hc_start, hc_A, 4)
        self.assertAlmostEqual(w.sum(), 3, places = 2)
        w = sparse_vaidya.generateWeight(hc_start, hc_A, 4)
        self.assertAlmostEqual(w.sum(), 4, places = 2)

        birk_start = np.array([0.5, 0.5, 0.5, 0.5])
        w = sparse_dikinls.generateWeight(birk_start, birk_A, 4)
        self.assertAlmostEqual(w.sum(), 1, places = 1)
        w = sparse_john.generateWeight(birk_start, birk_A, 4)
        self.assertAlmostEqual(w.sum(), 1.5, places = 2)
        w = sparse_vaidya.generateWeight(birk_start, birk_A, 4)
        self.assertAlmostEqual(w.sum(), 2, places = 2)


if __name__ == '__main__':
    unittest.main()