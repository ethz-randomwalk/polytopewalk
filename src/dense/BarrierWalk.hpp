#ifndef BARRIER_HPP
#define BARRIER_HPP

#include "RandomWalk.hpp"

class BarrierWalk : public RandomWalk{
    public:
        
        /**
         * @brief initialization of BarrierWalk class
         * @param r spread parameter
         * @param thin thin constant
         */
        BarrierWalk(double r, int thin = 1) : R(r), RandomWalk(thin){

        }

        /**
         * @brief generate weights when calculating Hessian matrix
         * @param x point in polytope to generate weight
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @return Vector
         */
        virtual VectorXd generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b);

        /**
         * @brief generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param init initial starting point
         * @param A polytope matrix
         * @param b polytope vector
         * @param burn number of initial steps to cut
         * @param seed seed for reproducibility
         * @return num_steps by d (dimension of x) matrix
         */
        MatrixXd generateCompleteWalk(const int num_steps, VectorXd& init, const MatrixXd& A, const VectorXd& b, int burn, int seed) override;

         /**
         * @brief set distribution constant
         * @param d (dimension)
         * @param n (number of constraints)
         */
        virtual void setDistTerm(int d, int n);
    
    protected:

        /**
         * @brief spread parameter
         */
        const double R;


        /**
         * @brief distribution constant
         */
        double DIST_TERM;

        /**
         * @brief generates a gaussian random vector with d dimension
         * @param d dimension
         * @param gen random number generator
         * @return Vector
         */
        VectorXd generateGaussianRV(int d, std::mt19937& gen);

        /**
         * @brief generates b - Ax (called slack) and 
         * makes global variable slack equal to it
         * @param x point
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @return vector
         */
        VectorXd generateSlack(const VectorXd& x, const MatrixXd& A, const VectorXd& b);

        /**
         * @brief calculates Mahalanobis distance weighted by Hessian matrix m
         * @param m Weighted Hessian Matrix
         * @param v vector to be measured
         * @return norm distance (double)
         */
        double localNorm(VectorXd v, const MatrixXd& m);

        /**
         * @brief generates Hessian of Log Barrier
         * @param x centered at x
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @return Matrix
         */
        MatrixXd generateHessian(const VectorXd& x, const MatrixXd& A, const VectorXd& b);

        /**
         * @brief generates a point drawn from a Multivariate Gaussian N(x, f(Hessian(x)))
         * @param x centered point in the polytope
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @param gen random number generator
         * @returns Vector
         */
        VectorXd generateSample(const VectorXd& x, const MatrixXd& A, const VectorXd& b, std::mt19937& gen);
};

#endif