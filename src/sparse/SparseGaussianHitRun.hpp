#ifndef CONSGENHITRUN_HPP
#define CONSGENHITRUN_HPP

#include "SparseRandomWalk.hpp"

class SparseGaussianHitAndRun : public SparseRandomWalk{
    public:
        /**
         * @brief initialization of Sparse Hit and Run class
         * @param r spread parameter
         * @param mu mean vector in Gaussian Distribution
         * @param cov covariance matrix in Gaussian Distribution
         * @param thin thin parameter
         * @param err error constant
         */
        SparseGaussianHitAndRun(double r, VectorXd mu, SparseMatrixXd cov, 
            int thin = 1, double err = 1e-6) : R(r), MU(mu), COV(cov), 
            SparseRandomWalk(thin, err) {
            
            chol.analyzePattern(COV);
            chol.factorize(COV);

        }

         /**
         * @brief generate values from the Hit and Run
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
        
    
    protected:
        /**
         * @brief spread parameter
         */
        const double R;

        /**
         * @brief mean of Gaussian Distribution
         */
        const VectorXd MU; 

        /**
         * @brief Covariance matrix of Gaussian Distribution
         */
        const SparseMatrixXd COV; 

        /**
         * @brief solve systems of equations with covariance matrix
         */
        SimplicialLLT<SparseMatrixXd> chol;

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
        double minF(VectorXd& v, VectorXd&x, double l, double u);

        /**
         * @brief runs binary search to find a suitable chord intersection with the polytope
         * @param direction (random direction variable)
         * @param x (starting point)
         * @param k k values >= 0 constraint
         * @return double 
         */
        double binarySearch(
            VectorXd direction, 
            VectorXd& x,
            int k);


};

#endif