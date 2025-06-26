import unittest
from polytopewalk import FacialReduction
import numpy as np 
import pandas as pd
from scipy.sparse import lil_matrix
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



class TestDenseWalk(unittest.TestCase):
    def test_run(self):
        dikin = DikinWalk(r = 0.5)
        vaidya = VaidyaWalk(r = 0.5)
        john = JohnWalk(r = 0.5)

        fr = FacialReduction()
        dc = DenseCenter()
        fr_res_simplex = fr.reduce(simplex_A, simplex_b, simplex_A.shape[1], False)
        fr_res_hc = fr.reduce(hc_A, hc_b, 4, False)
        fr_res_birk = fr.reduce(birk_A, birk_b, birk_A.shape[1], False)

        init_simplex = dc.getInitialPoint(fr_res_simplex.dense_A, fr_res_simplex.dense_b)
        init_hc = dc.getInitialPoint(fr_res_hc.dense_A, fr_res_hc.dense_b)
        init_birk = dc.getInitialPoint(fr_res_birk.dense_A, fr_res_birk.dense_b)

        res1 = dikin.generateCompleteWalk(100, init_simplex, fr_res_simplex.dense_A, fr_res_simplex.dense_b, 10, 10, seed = 5)
        self.assertTrue(res1.shape == (9, 2))
        res2 = vaidya.generateCompleteWalk(100, init_simplex, fr_res_simplex.dense_A, fr_res_simplex.dense_b, 10, 10, seed = 5)
        self.assertTrue(res2.shape == (9, 2))
        res3 = john.generateCompleteWalk(100, init_simplex, fr_res_simplex.dense_A, fr_res_simplex.dense_b, 10, 10, seed = 5)
        self.assertTrue(res3.shape == (9, 2))


        res1 = dikin.generateCompleteWalk(100, init_hc, fr_res_hc.dense_A, fr_res_hc.dense_b, 10, 10, seed = 5)
        self.assertTrue(res1.shape[0] == 9)
        res2 = vaidya.generateCompleteWalk(100, init_hc, fr_res_hc.dense_A, fr_res_hc.dense_b, 10, 10, seed = 5)
        self.assertTrue(res2.shape[0] == 9)
        res3 = john.generateCompleteWalk(100, init_hc, fr_res_hc.dense_A, fr_res_hc.dense_b, 10, 10, seed = 5)
        self.assertTrue(res3.shape[0] == 9)

        res1 = dikin.generateCompleteWalk(100, init_birk, fr_res_birk.dense_A, fr_res_birk.dense_b, 10, 10, seed = 5)
        self.assertTrue(res1.shape[0] == 9)
        res2 = vaidya.generateCompleteWalk(100, init_birk, fr_res_birk.dense_A, fr_res_birk.dense_b, 10, 10, seed = 5)
        self.assertTrue(res2.shape[0] == 9)
        res3 = john.generateCompleteWalk(100, init_birk, fr_res_birk.dense_A, fr_res_birk.dense_b, 10, 10, seed = 5)
        self.assertTrue(res3.shape[0] == 9)
 
if __name__ == '__main__':
    unittest.main()

