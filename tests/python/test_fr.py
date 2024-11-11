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

class TestFacialReduction(unittest.TestCase):

    def test_fr(self):
        fr = FacialReduction()
        simplex_dense = fr.reduce(simplex_A, simplex_b, 3, False)
        hc_dense = fr.reduce(hc_A, hc_b, 4, False)
        birk_dense = fr.reduce(birk_A, birk_b, 4, False)

        self.assertTrue(simplex_dense.dense_A.shape == (3,2))
        self.assertTrue(len(simplex_dense.dense_b) == 3)

        self.assertTrue(hc_dense.dense_A.shape == (4,2))
        self.assertTrue(len(hc_dense.dense_b) == 4)

        self.assertTrue(birk_dense.dense_A.shape == (4,1))
        self.assertTrue(len(birk_dense.dense_b) == 4)

        A = np.array([[1, 1, 0], [-1, -1, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])
        A = np.hstack((A, np.eye(6)))
        b = np.array([1, -1, 1, 1, 1, 1])

        fr_res = fr.reduce(A, b, 6, False)
        self.assertTrue(fr_res.sparse_A.shape == (5, 7))
        self.assertTrue(fr_res.dense_A.shape == (4, 2))


        A = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])
        A = np.hstack((A, np.eye(6)))
        b = np.array([1, 1, 0, 0, 0, 0])

        fr_res = fr.reduce(A, b, 6, False)
        self.assertTrue(fr_res.sparse_A.shape == (4, 5))
        self.assertTrue(fr_res.dense_A.shape == (2, 1))

        A = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])
        A = np.hstack((A, np.eye(4)))
        b = np.array([1, 0, 1, 0])

        fr_res = fr.reduce(A, b, 4, False)
        self.assertTrue(fr_res.sparse_A.shape == (4, 6))
        self.assertTrue(fr_res.dense_A.shape == (4, 2))


    
if __name__ == '__main__':
    unittest.main()

