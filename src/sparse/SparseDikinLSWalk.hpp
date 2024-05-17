#ifndef CONSDIKINLS_HPP
#define CONSDIKINLS_HPP

#include "SparseBarrierWalk.hpp"
#include "LeverageScore.hpp"

class SparseDikinLSWalk : public SparseBarrierWalk{

    public:
        /**
         * @brief constructor for weighted Dikin Walk class
         * @param r spread parameter
         * @param thin thin parameter
         * @param g_lim gradient descent norm limit
         * @param step_size size of gradient descent step
         * @param max_iter maximum number of iterations in gradient descent
         * @param err error constant
         */
        SparseDikinLSWalk(double r, int thin = 1, double g_lim = 0.01, double step_size = 0.1, int max_iter = 1000, double err = 1e-6) : G_LIM(g_lim), STEP_SIZE(step_size), MAX_ITER(max_iter), SparseBarrierWalk(r, thin, err) {}

        /**
         * @brief generate weight by solving convex optimization task
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
         * @return void
         */
        void setDistTerm(int d, int n) override;

        /**
         * @brief stops gradient descent if it reaches under this number
         */
        const double G_LIM;

        /**
         * @brief step size for gradient descent
         */
        const double STEP_SIZE;

        /**
         * @brief max number of iterations in gradient descent
         */
        const int MAX_ITER;

        /**
         * @brief saves current weight for iteration
         */
        VectorXd w_i = VectorXd::Zero(1) - VectorXd::Ones(1); 

};

#endif