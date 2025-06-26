import unittest
from polytopewalk import FacialReduction
import numpy as np 
import pandas as pd
from scipy.sparse import lil_matrix
from polytopewalk.sparse import *

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



class TestSparseWalk(unittest.TestCase):
    def test_run(self):
        dikin = SparseDikinWalk(r = 0.8)
        vaidya = SparseVaidyaWalk(r = 0.8)
        john = SparseJohnWalk(r = 0.8)

        fr = FacialReduction()
        sc = SparseCenter()
        fr_res_simplex = fr.reduce(simplex_A, simplex_b, simplex_A.shape[1], True)
        fr_res_hc = fr.reduce(hc_A, hc_b, 4, True)
        fr_res_birk = fr.reduce(birk_A, birk_b, birk_A.shape[1], True)

        init_simplex = sc.getInitialPoint(fr_res_simplex.sparse_A, fr_res_simplex.sparse_b, 3)
        init_hc = sc.getInitialPoint(fr_res_hc.sparse_A, fr_res_hc.sparse_b, 4)
        init_birk = sc.getInitialPoint(fr_res_birk.sparse_A, fr_res_birk.sparse_b, 4)

        res1 = dikin.generateCompleteWalk(100, init_simplex, fr_res_simplex.sparse_A, fr_res_simplex.sparse_b, 3, 10, 10, seed = 5)
        self.assertTrue(res1.shape == (9, 3))
        res2 = vaidya.generateCompleteWalk(100, init_simplex, fr_res_simplex.sparse_A, fr_res_simplex.sparse_b, 3, 10, 10, seed = 5)
        self.assertTrue(res2.shape == (9, 3))
        res3 = john.generateCompleteWalk(100, init_simplex, fr_res_simplex.sparse_A, fr_res_simplex.sparse_b, 3, 10, 10, seed = 5)
        self.assertTrue(res3.shape == (9, 3))


        res1 = dikin.generateCompleteWalk(100, init_hc, fr_res_hc.sparse_A, fr_res_hc.sparse_b, 4, 10, 10, seed = 5)
        self.assertTrue(res1.shape[0] == 9)
        res2 = vaidya.generateCompleteWalk(100, init_hc, fr_res_hc.sparse_A, fr_res_hc.sparse_b, 4, 10, 10, seed = 5)
        self.assertTrue(res2.shape[0] == 9)
        res3 = john.generateCompleteWalk(100, init_hc, fr_res_hc.sparse_A, fr_res_hc.sparse_b, 4, 10, 10, seed = 5)
        self.assertTrue(res3.shape[0] == 9)

        res1 = dikin.generateCompleteWalk(100, init_birk, fr_res_birk.sparse_A, fr_res_birk.sparse_b, 4, 10, 10, seed = 5)
        self.assertTrue(res1.shape[0] == 9)
        res2 = vaidya.generateCompleteWalk(100, init_birk, fr_res_birk.sparse_A, fr_res_birk.sparse_b, 4, 10, 10, seed = 5)
        self.assertTrue(res2.shape[0] == 9)
        res3 = john.generateCompleteWalk(100, init_birk, fr_res_birk.sparse_A, fr_res_birk.sparse_b, 4, 10, 10, seed = 5)
        self.assertTrue(res3.shape[0] == 9)
 
if __name__ == '__main__':
    unittest.main()
