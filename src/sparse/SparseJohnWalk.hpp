#ifndef CONSJOHN_HPP
#define CONSJOHN_HPP

#include "SparseBarrierWalk.hpp"

class SparseJohnWalk : public SparseBarrierWalk{

    public:
        /**
         * @brief initialization of Sparse John Walk class
         * @param r spread parameter
         * @param lambda hessian regularization term
         * @param dist_func convex function f such that x follows exp(-f(x)) log-concave dist
         * @param thin thin parameter
         * @param lim limit in l-infinity norm
         * @param max_iter maximum number of iterations in fixed iteration
         * @param err error constant
         */
        SparseJohnWalk(double r, double lambda, function<double(const VectorXd&)> dist_func,
            int thin = 1, double lim = 1e-5, int max_iter = 1000, double err = 1e-5) : LIM(lim),
            MAX_ITER(max_iter), SparseBarrierWalk(r, lambda, dist_func, thin, err) {}

        /**
         * @brief generate weight by solving fixed point iteration
         * @param x slack variable
         * @param A polytope constraint
         * @param k k values >= 0 constraint
         * @return SparseMatrixXd
         */
        SparseMatrixXd generateWeight(
            const VectorXd& x, 
            const SparseMatrixXd& A,
            int k
        ) override; 
    
    protected:

        /**
         * @brief set distribution constant
         * @param d polytope matrix 
         * @param n polytope vector
         */
        void setDistTerm(int d, int n) override;

        /**
         * @brief stops if it reaches under this number during fixed iteration
         */
        const double LIM;

        /**
         * @brief max number of iterations in fixed iteration
         */
        const int MAX_ITER;

        /**
         * @brief saves current weight for iteration
         */
        VectorXd w_i = VectorXd::Zero(1) - VectorXd::Ones(1); 



};

#endif