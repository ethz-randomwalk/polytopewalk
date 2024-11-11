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


class TestInitialization(unittest.TestCase):
    def test_init(self):
        sc = SparseCenter()
        dc = DenseCenter()
        
        x = sc.getInitialPoint(simplex_A, simplex_b, 3)
        self.assertAlmostEqual(max(abs(x - np.array([1/3, 1/3, 1/3]))), 0)
        
        x = sc.getInitialPoint(hc_A, hc_b, 4)
        self.assertAlmostEqual(max(abs(x - np.array([0, 0, 1, 1, 1, 1]))), 0)

        x = sc.getInitialPoint(birk_A, birk_b, 4)
        self.assertAlmostEqual(max(abs(x - np.array([0.5, 0.5, 0.5, 0.5]))), 0)

        A = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]])
        b = np.array([1,1,0,0])
        x = dc.getInitialPoint(A, b)
        self.assertAlmostEqual(max(abs(x - np.array([0.5, 0.5]))), 0)

        A = np.array([[-1, 0], [0, -1], [1, 1]])
        b = np.array([0, 0, 1])
        x = dc.getInitialPoint(A, b)
        self.assertAlmostEqual(max(abs(x - np.array([1/3, 1/3]))), 0)



if __name__ == '__main__':
    unittest.main()