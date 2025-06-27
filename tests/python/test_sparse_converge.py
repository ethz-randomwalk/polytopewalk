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

def generate_simplex(d):
    return np.array([1/d] * d), np.array([[1] * d]), np.array([1]), d


class TestWeights(unittest.TestCase):
    def test_sparse_converge(self):
        x1, A1, b1, k1 = generate_simplex(5)

        k = 5
        cov_simplex = np.full((k,k), (-1.0 / (k*k)) / (k + 1.0))
        for i in range(k):
            cov_simplex[i,i] = (1 - (1.0 / (k))) * (1.0/k) / (k + 1.0)

        nchains = 5
        niter = 100_000
        burn = 1000
        thin = 10
        
        seed = [1000, 2000, 3000, 4000, 5000]
        sparse_dikin = SparseDikinWalk(r = 0.5)
        sparse_vaidya = SparseVaidyaWalk(r = 0.7)
        sparse_john = SparseJohnWalk(r = 0.7)
        
        walks_all_dikin = np.empty([nchains, (niter-burn)//thin, k])
        walks_all_vaidya = np.empty([nchains, (niter-burn)//thin, k])
        walks_all_john = np.empty([nchains, (niter-burn)//thin, k])
        for c in range(nchains):
            walks_all_dikin[c, :] = sparse_dikin.generateCompleteWalk(niter, x1, A1, b1, k1, 
                                                                      burnin = burn, thin = thin, seed = seed[c])
            walks_all_vaidya[c, :] = sparse_vaidya.generateCompleteWalk(niter, x1, A1, b1, k1, 
                                                                        burnin = burn, thin = thin, seed = seed[c])
            walks_all_john[c, :] = sparse_john.generateCompleteWalk(niter, x1, A1, b1, k1, 
                                                                    burnin = burn, thin = thin, seed = seed[c])

        walks_dikin = arviz.convert_to_dataset(walks_all_dikin)
        max_rhat_dikin = arviz.rhat(walks_dikin).to_array().values.max().item()
        
        walks_vaidya = arviz.convert_to_dataset(walks_all_vaidya)
        max_rhat_vaidya = arviz.rhat(walks_vaidya).to_array().values.max().item()
        
        walks_john = arviz.convert_to_dataset(walks_all_john)
        max_rhat_john = arviz.rhat(walks_john).to_array().values.max().item()
        
        max_rhat_dikin, max_rhat_vaidya, max_rhat_john

        self.assertTrue(max_rhat_dikin < 1.01)
        self.assertTrue(max_rhat_vaidya < 1.01)
        self.assertTrue(max_rhat_john < 1.01)

        walks_all_dikin_np = np.concatenate([walks_all_dikin[i,:,:] for i in range(nchains)])
        dikin_cov = np.cov(walks_all_dikin_np, rowvar = False)
        
        walks_all_vaidya_np = np.concatenate([walks_all_vaidya[i,:,:] for i in range(nchains)])
        vaidya_cov = np.cov(walks_all_vaidya_np, rowvar = False)
        
        walks_all_john_np = np.concatenate([walks_all_john[i,:,:] for i in range(nchains)])
        john_cov = np.cov(walks_all_john_np, rowvar = False)
        
        dikin_cov_diff = np.linalg.norm(dikin_cov - cov_simplex, ord='fro')
        vaidya_cov_diff = np.linalg.norm(vaidya_cov - cov_simplex, ord='fro')
        john_cov_diff = np.linalg.norm(john_cov - cov_simplex, ord='fro')
        
        dikin_mean_diff = np.linalg.norm(walks_all_dikin_np.mean(axis = 0) - [0.2] * 5)
        vaidya_mean_diff = np.linalg.norm(walks_all_vaidya_np.mean(axis = 0) - [0.2] * 5)
        john_mean_diff = np.linalg.norm(walks_all_john_np.mean(axis = 0) - [0.2] * 5)

        self.assertTrue(dikin_cov_diff < 5e-2)
        self.assertTrue(vaidya_cov_diff < 5e-2)
        self.assertTrue(john_cov_diff < 5e-2)

        self.assertTrue(dikin_mean_diff < 5e-2)
        self.assertTrue(vaidya_mean_diff < 5e-2)
        self.assertTrue(john_mean_diff < 5e-2)
        

if __name__ == '__main__':
    unittest.main()