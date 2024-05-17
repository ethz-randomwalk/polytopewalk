#ifndef CONSTRAINTWALK_HPP
#define CONSTRAINTWALK_HPP
#include "Common.hpp"

class SparseRandomWalk{

    public:
        /**
         * @brief Random Walk Superclass implementation
         * @param thin thin constant
         * @param err error constant
         */
        SparseRandomWalk(int thin = 1, double err = 1e-6) : ERR(err), THIN(thin) {}
    
        /**
         * @brief Generate values from the RandomWalk
         * @param num_steps number of steps wanted to take
         * @param init initial starting point
         * @param A polytope matrix 
         * @param b polytope vector
         * @param k k values >= 0 constraint
         * @return Matrix
         */
        virtual MatrixXd generateCompleteWalk(
            const int num_steps, 
            const VectorXd& init, 
            const SparseMatrixXd& A, 
            const VectorXd& b, 
            int k, 
            int burn
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
         * @return vector 
         */
        VectorXd generateGaussianRV(const int d);

        /**
         * @brief error constant 
         */
        const double ERR; 

        /**
         * @brief include every ___ sample
         */
        const int THIN; 
};

#endif