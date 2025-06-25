#ifndef RANDOMWALK_HPP
#define RANDOMWALK_HPP
#include "Common.hpp"

class RandomWalk{

    public:
    
        /**
         * @brief initialization of Random Walk super class
         * @param thin thin constant
         */
        RandomWalk(int thin = 1) : THIN(thin){}

        /**
         * @brief generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param init initial starting point
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @param burn number of initial steps to cut
         * @param seed seed for reproducibility
         * @return num_steps by d (dimension of x) matrix
         */
        virtual MatrixXd generateCompleteWalk(const int num_steps, VectorXd& init, const MatrixXd& A, const VectorXd& b, int burn, int seed);

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
         * @brief only include every __ sample
         */
        const int THIN;


        /**
         * @brief initialize random number generator
         * @param seed seed number for reproducible results
         * @return mt19937 random number generator
         */
        std::mt19937 initializeRNG(int seed);

};

#endif