#ifndef RANDOMWALK_HPP
#define RANDOMWALK_HPP
#include "Common.hpp"

class RandomWalk{

    public:
    
        /**
         * @brief initialization of Random Walk super class
         */
        RandomWalk(){}

        /**
         * @brief generate values from the walk
         * @param niter number of steps wanted to take
         * @param init initial starting point
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @param burnin number of initial steps to cut
         * @param thin thinning parameter
         * @param seed seed for reproducibility
         * @return (niter - burnin)//thin by d (dimension of x) matrix
         */
        virtual MatrixXd generateCompleteWalk(const int niter, VectorXd& init, const MatrixXd& A, const VectorXd& b, 
            int burnin, int thin, int seed);

    protected: 

        /**
         * @brief checks Az <= b
         * @param z vector
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @return bool (inside polytope or not)
         */
        bool inPolytope(const VectorXd& z, const MatrixXd& A, const VectorXd& b);

        /**
         * @brief returns normalized Gaussian vector of dimension d
         * @param d
         * @param gen random number generator
         * @return vector (normalized vector)
         */
        VectorXd generateGaussianRVNorm(const int d, std::mt19937& gen);

        /**
         * @brief prints unique identifier of the walk
         */
        virtual void printType();

        /**
         * @brief initialize random number generator
         * @param seed seed number for reproducible results
         * @return mt19937 random number generator
         */
        std::mt19937 initializeRNG(int seed);

};

#endif