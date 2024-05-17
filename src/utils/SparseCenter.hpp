#ifndef SPARSE_CENTER_HPP
#define SPARSE_CENTER_HPP

#include "Common.hpp"
#include "SparseLP.hpp"

class SparseCenter {
    public:
        /**
         * @brief Center Algorithm
         * @param max_iter maximum iterations during linear program
         * @param tol tolerance term
         * @param s_max error term
         * @param err_lp error for LP problem
         * @return res
         */
        SparseCenter(int max_iter = 3000, double tol = 1e-8, double s_max = 100.0, double err_lp = 1e-8) : MAX_ITER(max_iter), TOL(tol), S_MAX(s_max), ERR_LP(err_lp) {

        }

        /**
         * @brief Finds analytical center Ax = b, x >=_k 0 
         * @param A polytope matrix (Ax = b)
         * @param b polytope vector (Ax = b)
         * @param k k values >= 0 constraint
         * @return VectorXd 
         */
        VectorXd getInitialPoint(SparseMatrixXd A, VectorXd b, int k);

    protected:
        /**
         * @brief max iterations on linear program
         */
        const int MAX_ITER;

        /**
         * @brief tolerance parameter
         */
        const double TOL;

        /**
         * @brief tolerance parameter
         */
        const double S_MAX;

        /**
         * @brief error in linear program tolerance
         */
        const double ERR_LP;
};

#endif