
#ifndef DIKINWALK_HPP
#define DIKINWALK_HPP

#include "BarrierWalk.hpp"

class DikinWalk: public BarrierWalk{

    public:
        /**
         * @brief initialization of Dikin Walk class
         * @param r spread parameter
         * @param lambda hessian regularization term
         * @param dist_func convex function f such that x follows exp(-f(x)) log-concave dist
         * @param thin thin parameter
         */
        DikinWalk(double r, double lambda = 0, function<double(const VectorXd&)> dist_func = [](const VectorXd& x) { return 1.0; }, 
            int thin = 1) : BarrierWalk(r, lambda, dist_func, thin){}
        
        /**
         * @brief print dikin
         */
        void printType() override;

        /**
         * @brief returns weight for DikinWalk (Identity Matrix)
         * @param x point
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         */
        void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd&b) override;

    protected:

        /**
         * @brief set distribution constant
         * @param d (dimension)
         * @param n (number of constraints)
         */
        void setDistTerm(int d, int n) override;

};


#endif