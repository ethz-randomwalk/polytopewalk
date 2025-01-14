
#ifndef VAIDYAWALK_HPP
#define VAIDYAWALK_HPP

#include "BarrierWalk.hpp"

class VaidyaWalk: public BarrierWalk{

    public:
        /**
         * @brief initialization of Sparse Vaidya Walk class
         * @param r spread parameter
         * @param lambda hessian regularization term
         * @param dist_type distribution type {uniform, normal, log-concave}
         * @param thin thin constant
         */
        VaidyaWalk(double r, double lambda = 0, string dist_type = "uniform", 
            int thin = 1) : BarrierWalk(r, lambda, dist_type, thin){}  

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