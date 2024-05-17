
#ifndef DIKINLSWALK_HPP
#define DIKINLSWALK_HPP

#include "BarrierWalk.hpp"

class DikinLSWalk: public BarrierWalk{

    public:
        /**
         * @brief constructor for weighted Dikin Walk class
         * @param r spread parameter
         * @param thin thin constant
         * @param g_lim gradient descent norm limit
         * @param step_size size of gradient descent step
         * @param max_iter maximum number of iterations in gradient descent
         */
        DikinLSWalk(double r, int thin = 1, double g_lim = 0.01, double step_size = 0.1, int max_iter = 1000) : STEPSIZE(step_size), MAXITER(max_iter), GRADLIM(g_lim), BarrierWalk(r, thin){
            
        }

        /**
         * @brief print dikinls
         * @return void
         */
        void printType() override;

        /**
         * @brief generate weights when calculating Hessian matrix
         * @param x point in polytope to generate DikinLS weight
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @return void (update global variable weights)
         */
        void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b) override;
    
    protected:
        /**
         * @brief step size for gradient descent
         */
        const double STEPSIZE;

        /**
         * @brief max number of iterations in gradient descent
         */
        const int MAXITER;

        /**
         * @brief stops gradient descent if it reaches under this number
         */
        const double GRADLIM;

        /**
         * @brief saves current weight for iteration
         */
        VectorXd w_i = VectorXd::Zero(1) - VectorXd::Ones(1); 

        /**
         * @brief set Distribution Constant
         * @param d (dimension)
         * @param n (number of constraints)
         * @return void
         */
        void setDistTerm(int d, int n) override;
};

#endif