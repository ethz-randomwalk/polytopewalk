#ifndef CONSBARRIERWALK_HPP
#define CONSBARRIERWALK_HPP
#include "SparseRandomWalk.hpp"

class SparseBarrierWalk : public SparseRandomWalk{

    public:
        /**
         * @brief SparseBarrierWalk class
         * @param r spread parameter
         * @param thin thin parameter
         * @param err error term parameter
         */
        SparseBarrierWalk(double r, int thin = 1, double err = 1e-6) : R(r), SparseRandomWalk(thin, err) {}

        /**
         * @brief generate weight for slack inverse
         * @param x slack variable
         * @param A polytope constraint
         * @param k k values >= 0 constraint
         * @return SparseMatrixXd
         */
        virtual SparseMatrixXd generateWeight(
            const VectorXd& x, 
            const SparseMatrixXd& A,
            int k
        );

        /**
         * @brief Generate values from the SparseBarrierWalk
         * @param num_steps number of steps wanted to take
         * @param init initial starting point
         * @param A polytope matrix 
         * @param b polytope vector
         * @param k k values >= 0 constraint
         * @param burn number of initial steps to cut
         * @return Matrix
         */
        MatrixXd generateCompleteWalk(
            const int num_steps, 
            const VectorXd& init, 
            const SparseMatrixXd& A, 
            const VectorXd& b, 
            int k,
            int burn
        ) override;
        
        /**
         * @brief Distribution constant
         * @param d polytope matrix 
         * @param n polytope vector
         * @return void
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
         * @return VectorXd
         */
        VectorXd generateSample(
            const VectorXd& x, 
            const SparseMatrixXd& A, 
            int k
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




