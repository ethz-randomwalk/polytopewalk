#ifndef CONSBALLWALK_HPP
#define CONSBALLWALK_HPP

#include "SparseRandomWalk.hpp"

class SparseBallWalk : public SparseRandomWalk{
    public:
        /**
         * @brief initialization of Sparse Ball Walk class
         * @param r spread parameter
         */
        SparseBallWalk(double r) : R(r), SparseRandomWalk(0.0){}

         /**
         * @brief generate values from the Ball walk (constrained)
         * @param niter number of steps wanted to take
         * @param init initial starting point
         * @param A polytope matrix 
         * @param b polytope vector
         * @param k k values >= 0 constraint
         * @param burnin number of initial steps to cut
         * @param thin thinning parameter
         * @param seed seed for reproducibility
         * @return (niter - burnin)//thin by d (dimension of x) matrix
         */
        MatrixXd generateCompleteWalk(
            const int niter, 
            const VectorXd& init, 
            const SparseMatrixXd& A, 
            const VectorXd& b, 
            int k, 
            int burnin,
            int thin,
            int seed
            ) override;
    
    protected:
        /**
         * @brief spread parameter
         */
        const double R;
};
#endif 