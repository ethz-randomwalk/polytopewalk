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
         * @param x initial starting point
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @param burn number of initial steps to cut
         * @return num_steps by d (dimension of x) matrix
         */
        virtual MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b, int burn);

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
         * @return vector (normalized vector)
         */
        VectorXd generateGaussianRVNorm(const int d);

        /**
         * @brief prints unique identifier of the walk
         */
        virtual void printType();

        /**
         * @brief only include every __ sample
         */
        const int THIN;
};

#endif