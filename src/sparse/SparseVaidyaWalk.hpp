#ifndef CONSVAIDYA_HPP
#define CONSVAIDYA_HPP

#include "SparseBarrierWalk.hpp"

class SparseVaidyaWalk : public SparseBarrierWalk{

    public:
        /**
         * @brief constructor for Vaidya Walk class
         * @param r spread parameter
         * @param lambda hessian regularization term
         * @param dist_func convex function f such that x follows exp(-f(x)) log-concave dist
         * @param thin thin parameter
         * @param err error constant
         */
        SparseVaidyaWalk(double r,double lambda, function<double(const VectorXd&)> dist_func,
            int thin = 1, double err = 1e-6) : SparseBarrierWalk(r, lambda, dist_func, thin, err) {}

        /**
         * @brief generate weight (leverage score calculation)
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
         * @brief Distribution constant
         * @param d polytope matrix 
         * @param n polytope vector
         */
        void setDistTerm(int d, int n) override;

};

#endif