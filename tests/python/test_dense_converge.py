import unittest
import numpy as np
import pandas as pd 
import arviz
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix, lil_matrix, csr_array

#Import Sparse Barrier Random Walks
from polytopewalk.sparse import SparseDikinWalk, SparseVaidyaWalk, SparseJohnWalk

#Import Dense Barrier Random Walks
from polytopewalk.dense import DikinWalk, VaidyaWalk, JohnWalk

class TestWeights(unittest.TestCase):
    def test_dense_converge(self):
        x2 = np.array([0] * 3)
        A2 = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])
        b2 = np.array([1] * 6)
        cov_hc = 1/3 * np.eye(3)

        nchains = 5
        niter = 300_000
        burn = 1_000
        thin = 10
        
        seed = [1000, 2000, 3000, 4000, 5000]
        dikin = DikinWalk(r = 0.8)
        vaidya = VaidyaWalk(r = 0.9)
        john = JohnWalk(r = 0.9)
        
        walks_all_dikin = np.empty([nchains, (niter-burn)//thin, 3])
        walks_all_vaidya = np.empty([nchains, (niter-burn)//thin, 3])
        walks_all_john = np.empty([nchains, (niter-burn)//thin, 3])
        for c in range(nchains):
            walks_all_dikin[c, :] = dikin.generateCompleteWalk(niter, x2, A2, b2, 
                                                               burnin = burn, thin = thin, seed = seed[c])
            walks_all_vaidya[c, :] = vaidya.generateCompleteWalk(niter, x2, A2, b2, 
                                                                 burnin = burn, thin = thin, seed = seed[c])
            walks_all_john[c, :] = john.generateCompleteWalk(niter, x2, A2, b2, 
                                                             burnin = burn, thin = thin, seed = seed[c])

        walks_dikin = arviz.convert_to_dataset(walks_all_dikin)
        max_rhat_dikin = arviz.rhat(walks_dikin).to_array().values.max().item()
        
        walks_vaidya = arviz.convert_to_dataset(walks_all_vaidya)
        max_rhat_vaidya = arviz.rhat(walks_vaidya).to_array().values.max().item()
        
        walks_john = arviz.convert_to_dataset(walks_all_john)
        max_rhat_john = arviz.rhat(walks_john).to_array().values.max().item()

        self.assertTrue(max_rhat_dikin < 1.01)
        self.assertTrue(max_rhat_vaidya < 1.01)
        self.assertTrue(max_rhat_john < 1.01)

        walks_all_dikin_np = np.concatenate([walks_all_dikin[i,:,:] for i in range(nchains)])
        dikin_cov = np.cov(walks_all_dikin_np, rowvar = False)
        
        walks_all_vaidya_np = np.concatenate([walks_all_vaidya[i,:,:] for i in range(nchains)])
        vaidya_cov = np.cov(walks_all_vaidya_np, rowvar = False)
        
        walks_all_john_np = np.concatenate([walks_all_john[i,:,:] for i in range(nchains)])
        john_cov = np.cov(walks_all_john_np, rowvar = False)
        
        dikin_cov_diff = np.linalg.norm(dikin_cov - cov_hc, ord='fro')
        vaidya_cov_diff = np.linalg.norm(vaidya_cov - cov_hc, ord='fro')
        john_cov_diff = np.linalg.norm(john_cov - cov_hc, ord='fro')
        
        dikin_mean_diff = np.linalg.norm(walks_all_dikin_np.mean(axis = 0) - [0.0] * 3)
        vaidya_mean_diff = np.linalg.norm(walks_all_vaidya_np.mean(axis = 0) - [0.0] * 3)
        john_mean_diff = np.linalg.norm(walks_all_john_np.mean(axis = 0) - [0.0] * 3)

        self.assertTrue(dikin_cov_diff < 5e-2)
        self.assertTrue(vaidya_cov_diff < 5e-2)
        self.assertTrue(john_cov_diff < 5e-2)

        self.assertTrue(dikin_mean_diff < 5e-2)
        self.assertTrue(vaidya_mean_diff < 5e-2)
        self.assertTrue(john_mean_diff < 5e-2)
                        

if __name__ == '__main__':
    unittest.main()