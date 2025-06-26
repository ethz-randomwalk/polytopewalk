
#ifndef HITRUN_HPP
#define HITRUN_HPP

#include "RandomWalk.hpp"

class HitAndRun: public RandomWalk{

    public:
        /**
         * @brief initialization of Hit and Run class
         * @param r spread hyperparamter
         * @param err error hyperparameter
         */
        HitAndRun(double r, double err = 1e-6) : ERR(err), R(r), RandomWalk() {

        }

        /**
         * @brief Generate values from the walk
         * @param niter number of steps wanted to take
         * @param init initial starting point
         * @param A polytope matrix
         * @param b polytope matrix
         * @param burnin number of steps to burn
         * @param thin thinning parameter
         * @param seed seed for reproducibility
         * @return (niter - burnin)//thin by d (dimension of x) matrix
         */
        MatrixXd generateCompleteWalk(const int niter, VectorXd& init, const MatrixXd& A, const VectorXd& b, int burnin, int thin, int seed) override;

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