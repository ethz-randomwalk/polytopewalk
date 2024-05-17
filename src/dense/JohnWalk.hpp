
#ifndef JOHNWALK_HPP
#define JOHNWALK_HPP

#include "BarrierWalk.hpp"

class JohnWalk: public BarrierWalk{

    public:
         /**
         * @brief constructor for John Walk class
         * @param r spread parameter
         * @param thin thin constant
         * @param lim norm limit for fixed point iteration
         * @param max_iter maximum number of iterations in gradient descent
         */
        JohnWalk(double r, int thin = 1, double lim = 1e-5, int max_iter = 1000) : MAXITER(max_iter), LIM(lim), BarrierWalk(r, thin){

        }

        /**
         * @brief print john walk
         * @return void
         */
        void printType() override;

        /**
         * @brief generates John weight by solving convex optimization problem
         * @param x point in polytope to generate DikinLS weight
         * @param A polytope matrix
         * @param b polytope matrix
         * @return void (update global variable weights)
         */
        void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b) override;

    
    protected:

        /**
         * @brief max number of iterations in fixed point iteration
         */
        const double MAXITER;

        /**
         * @brief stops if it reaches under this number in fixed iteration
         */
        const double LIM;

        /**
         * @brief saves current weight for iteration
         */
        VectorXd w_i = VectorXd::Zero(1) - VectorXd::Ones(1); 

        /**
         * @brief set Dist Term for John Walk
         * @param d (dimension)
         * @param n (number of constraints)
         * @return void
         */
        void setDistTerm(int d, int n) override;

};

#endif