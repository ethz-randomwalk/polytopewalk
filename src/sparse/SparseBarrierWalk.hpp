#ifndef CONSBARRIERWALK_HPP
#define CONSBARRIERWALK_HPP
#include "SparseRandomWalk.hpp"

class SparseBarrierWalk : public SparseRandomWalk{

    public:
        /**
         * @brief initialization of Sparse Barrier Walk class
         * @param r spread parameter
         * @param err error term parameter
         */
        SparseBarrierWalk(double r, double err = 1e-6) : R(r), SparseRandomWalk(err) {}

        /**
         * @brief generate weight for slack inverse
         * @param x slack variable
         * @param A polytope constraint
         * @param k k values >= 0 constraint
         * @return Vector
         */
        virtual VectorXd generateWeight(
            const VectorXd& x, 
            const SparseMatrixXd& A,
            int k
        );

        /**
         * @brief generate values from the SparseBarrierWalk
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
        
        /**
         * @brief set distribution constant
         * @param d polytope matrix 
         * @param n polytope vector
         */
        virtual void setDistTerm(int d, int n); 
    
    protected:
        /**
         * @brief distribution constant
         */
        double DIST_TERM; 

        /**
         * @brief spread parameter
         */
        double R; 

        /**
         * @brief inverse solver
         */
        SparseLU<SparseMatrixXd> A_solver;

        /**
         * @brief generate slack inverse (1/x)
         * @param x vector value
         * @param k k values >= 0 constraint
         * @return SparseMatrixXd
         */
        SparseMatrixXd generateSlackInverse(
            const VectorXd& x, 
            int k
        );

        /**
         * @brief generate sample from distribution
         * @param x vector value
         * @param A polytope matrix (Ax = b)
         * @param k values >= 0 constraint
         * @param gen random number generator
         * @return VectorXd
         */
        VectorXd generateSample(
            const VectorXd& x, 
            const SparseMatrixXd& A, 
            int k,
            std::mt19937& gen
        ); 
        
        /**
         * @brief generate density term
         * @param x center of distribution
         * @param z value from distribution
         * @param A polytope matrix (Ax = b)
         * @param k values >= 0 constraint
         * @return double
         */
        double generateProposalDensity(
            const VectorXd& x, 
            const VectorXd& z, 
            const SparseMatrixXd& A, 
            int k
        );
};

#endif




