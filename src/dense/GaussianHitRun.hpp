
#ifndef GENHITRUN_HPP
#define GENHITRUN_HPP

#include "RandomWalk.hpp"

class GaussianHitAndRun: public RandomWalk{

    public:
        /**
         * @brief initialization of General Hit and Run class
         * @param r spread hyperparameter
         * @param mu mean vector in Gaussian Distribution
         * @param cov covariance matrix in Gaussian Distribution
         * @param dist_func distribution
         * @param thin thin parameter (record every ith value)
         * @param err error hyperparameter
         */
        GaussianHitAndRun(double r, VectorXd mu, MatrixXd cov, int thin = 1,
        double err = 1e-6) : R(r), MU(mu), COV(cov), ERR(err), RandomWalk(thin) {
            cov_inv = cov.inverse();
        }

        /**
         * @brief Generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @param A polytope matrix
         * @param b polytope matrix
         * @param burn number of steps to burn
         * @return num_steps by d (dimension of x) matrix
         */
        MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b, int burn) override;

         /**
         * @brief print general type 
         */
        void printType() override;
    
    protected:
        /**
         * @brief relative error of the binary search operation
         */
        const double ERR;

        /**
         * @brief mean of Gaussian Distribution
         */
        const VectorXd MU; 

        /**
         * @brief Covariance matrix of Gaussian Distribution
         */
         const MatrixXd COV; 


         /**
         * @brief Inverse Covariance matrix of Gaussian Distribution
         */
         MatrixXd cov_inv; 

        /**
         * @brief initial starting value
         */
        const double R;

        /**
         * @brief get distance between vectors x and y
         * @param x
         * @param y
         * @return double
         */
        double distance(VectorXd& x, VectorXd&y);

        /**
         * @brief get the gaussian log PDF 
         * @param x vector
         * @return double
         */
        double gaussianLogPDF(VectorXd x);
        
        
        /**
         * @brief get min f(x) in direction v through point x
         * @param v direction vector
         * @param x starting point
         * @param l lower bound
         * @param u upper bound
         * @return double
         */
        double minF(const VectorXd& v, const VectorXd&x, double l, double u);

        /**
         * @brief runs binary search to find a suitable chord intersection with the polytope
         * @param direction (random direction variable)
         * @param x (starting point)
         * @param A polytope matrix
         * @param b polytope vector
         * @return double 
         */
        double binarySearch(VectorXd direction, VectorXd& x, const MatrixXd& A, const VectorXd& b);

};

#endif