#ifndef CPF_HPP
#define CPF_HPP

#include "Common.hpp"


class DenseCenter {
    public:

        /**
         * @brief Center Algorithm
         * @param max_iter maximum iterations during linear program
         * @param tol tolerance term
         * @param s_max error term
         * @param start starting constant in linear program
         * @return res
         */
        DenseCenter(int max_iter = 3000, double tol = 1e-8, double s_max = 100, int start = 100000) : MAX_ITER(max_iter), TOL(tol), S_MAX(s_max), START(start) {};

        /**
         * @brief Finds analytical center Ax <= b
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @return VectorXd 
         */
        VectorXd getInitialPoint(MatrixXd A, VectorXd b);

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
         * @brief starting parameter
         */
        const int START;

};

#endif