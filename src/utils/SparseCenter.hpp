#ifndef SPARSE_CENTER_HPP
#define SPARSE_CENTER_HPP

#include "Common.hpp"
#include "SparseLP.hpp"

class SparseCenter {
    public:
        /**
         * @brief Center Algorithm for sparse input
         * @return res
         */
        SparseCenter(){};

        /**
         * @brief Finds analytical center Ax = b, x >=_k 0 
         * @param A polytope matrix (Ax = b)
         * @param b polytope vector (Ax = b)
         * @param k k values >= 0 constraint
         * @return VectorXd 
         */
        VectorXd getInitialPoint(SparseMatrixXd& A, VectorXd& b, int k);

};

#endif