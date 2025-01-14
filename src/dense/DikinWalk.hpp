
#ifndef DIKINWALK_HPP
#define DIKINWALK_HPP

#include "BarrierWalk.hpp"

class DikinWalk: public BarrierWalk{

    public:
        /**
         * @brief initialization of Dikin Walk class
         * @param r spread parameter
         * @param lambda hessian regularization term
         * @param dist_type distribution type {uniform, normal, log-concave}
         * @param thin thin parameter
         */
        DikinWalk(double r, double lambda = 0, string dist_type = "uniform", 
            int thin = 1) : BarrierWalk(r, lambda, dist_type, thin){}
        
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