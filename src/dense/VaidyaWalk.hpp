
#ifndef VAIDYAWALK_HPP
#define VAIDYAWALK_HPP

#include "BarrierWalk.hpp"

class VaidyaWalk: public BarrierWalk{

    public:
        /**
         * @brief initialization of Sparse Vaidya Walk class
         * @param r spread parameter
         * @param lambda hessian regularization term
         * @param dist_func convex function f such that x follows exp(-f(x)) log-concave dist
         * @param thin thin constant
         */
        VaidyaWalk(double r, double lambda = 0, function<double(const VectorXd&)> dist_func = [](const VectorXd& x) { return 1.0; }, 
            int thin = 1) : BarrierWalk(r, lambda, dist_func, thin){}  

        /**
         * @brief print general type 
         */
        void printType() override;

        /**
         * @brief returns weight for Vaidya Walk (leverage score calculation)
         * @param x center vector
         * @param A polytope matrix
         * @param b polytope vector
         */
        void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b) override;
    
    protected:

        /**
         * @brief set distribution constant
         * @param d (dimension)
         * @param n (number of constraints)
         */
        void setDistTerm(int d, int n) override;

};

#endif