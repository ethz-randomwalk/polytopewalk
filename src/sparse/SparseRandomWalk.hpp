#ifndef CONSTRAINTWALK_HPP
#define CONSTRAINTWALK_HPP
#include "Common.hpp"

class SparseRandomWalk{

    public:
        /**
         * @brief initialization of Sparse Random Walk class
         * @param err error constant
         */
        SparseRandomWalk(double err = 1e-6) : ERR(err){}
    
        /**
         * @brief Generate values from the RandomWalk
         * @param num_steps number of steps wanted to take
         * @param init initial starting point
         * @param A polytope matrix 
         * @param b polytope vector
         * @param k k values >= 0 constraint
         * @param burnin number of steps to burn
         * @param thin thinning parameter
         * @param seed seed for reproducibility
         * @return (niter - burnin)//thin by d (dimension of x) matrix
         */
        virtual MatrixXd generateCompleteWalk(
            const int num_steps, 
            const VectorXd& init, 
            const SparseMatrixXd& A, 
            const VectorXd& b, 
            int k, 
            int burnin,
            int thin,
            int seed
            );
        
    protected:
        /**
         * @brief check if value is in polytope
         * @param z proposal vector (assuming sampled from Ax = 0)
         * @param k k values >= 0 constraint
         * @return Matrix
         */
        bool inPolytope(const VectorXd& z, int k);

        /**
         * @brief returns Gaussian vector of dimension d
         * @param d
         * @param gen random number generator
         * @return vector 
         */
        VectorXd generateGaussianRV(const int d, std::mt19937& gen);


        /**
         * @brief initialize random number generator
         * @param seed seed number for reproducible results
         * @return mt19937 random number generator
         */
        std::mt19937 initializeRNG(int seed);

        /**
         * @brief error constant 
         */
        const double ERR; 
};

#endif