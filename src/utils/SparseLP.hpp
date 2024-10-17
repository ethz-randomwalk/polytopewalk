#ifndef SPARSE_LP_HPP
#define SPARSE_LP_HPP

#include "Common.hpp"

class SparseLP{
    public:
        /**
         * @brief Sparse Linear Programming Solver
         */
        SparseLP(){};

        /**
         * @brief Finds analytical center Ax = b, x >=_k 0 
         * @param A constraint matrix
         * @param b constriant vector
         * @param c objective vector 
         * @param row_rel relationship for A and b
         * @param col_cons constraint for columns
         * @param col_rel relation for columns
         * @return VectorXd 
         */
        VectorXd findOptimalVector(SparseMatrixXd& A, VectorXd& b, VectorXd& c, VectorXd& row_rel, 
                                   VectorXd& col_cons, VectorXd& col_rel); 

};

#endif