#ifndef CONSHITRUN_HPP
#define CONSHITRUN_HPP

#include "SparseRandomWalk.hpp"

class SparseHitAndRun : public SparseRandomWalk{
    public:
        /**
         * @brief initialization of Sparse Hit and Run class
         * @param r spread parameter
         * @param err error constant
         */
        SparseHitAndRun(double r, double err = 1e-6) : R(r), SparseRandomWalk(err) {}

         /**
         * @brief generate values from the Hit and Run
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

        /**
         * @brief runs binary search to find a suitable chord intersection with the polytope
         * @param direction (random direction variable)
         * @param x (starting point)
         * @param k k values >= 0 constraint
         * @return double 
         */
        double binarySearch(
            VectorXd direction, 
            VectorXd& x,
            int k);
};

#endif