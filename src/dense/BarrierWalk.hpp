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
         * @brief weights generated from generateWeights function
         */
        DiagonalMatrix<double, Dynamic> weights{};

        /**
         * @brief generate weights when calculating Hessian matrix
         * @param x point in polytope to generate weight
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         */
        virtual void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b);

        /**
         * @brief generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @param A polytope matrix
         * @param b polytope vector
         * @param burn number of initial steps to cut
         * @return num_steps by d (dimension of x) matrix
         */
        MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b, int burn) override;

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
         * @brief represents global variable b - Ax
         */
        VectorXd slack{}; 

        /**
         * @brief hessian matrix from global variable from generateHessian
         */
        MatrixXd hess{};

        /**
         * @brief new proposal point generated from generateSample function
         */
        VectorXd prop{};

        /**
         * @brief generates a gaussian random vector with d dimension
         * @param d dimension
         * @return Vector
         */
        VectorXd generateGaussianRV(int d);

        /**
         * @brief generates b - Ax (called slack) and 
         * makes global variable slack equal to it
         * @param x point
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         */
        void generateSlack(const VectorXd& x, const MatrixXd& A, const VectorXd& b);

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
         */
        void generateHessian(const VectorXd& x, const MatrixXd& A, const VectorXd& b);

        /**
         * @brief generates a point drawn from a Multivariate Gaussian N(x, f(Hessian(x)))
         * @param x centered point in the polytope
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         */
        void generateSample(const VectorXd& x, const MatrixXd& A, const VectorXd& b);
};

#endif