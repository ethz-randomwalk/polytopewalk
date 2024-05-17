#ifndef RANDOMWALK_HPP
#define RANDOMWALK_HPP
#include "Common.hpp"

class RandomWalk{

    public:
    
        /**
         * @brief Random Walk Superclass implementation
         * @param thin thin constant
         */
        RandomWalk(int thin = 1) : THIN(thin){}

        /**
         * @brief Generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @param burn number of initial steps to cut
         * @return Matrix
         */
        virtual MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b, int burn);

    protected: 

        /**
         * @brief checks Az <= b
         * @param z vector
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @return void
         */
        bool inPolytope(const VectorXd& z, const MatrixXd& A, const VectorXd& b);

        /**
         * @brief returns normalized Gaussian vector of dimension d
         * @param d
         * @return vector 
         */
        VectorXd generateGaussianRVNorm(const int d);

        /**
         * @brief prints unique identifier of the walk
         * @return void
         */
        virtual void printType();

        /**
         * @brief only include every __ sample
         */
        const int THIN;
};

#endif