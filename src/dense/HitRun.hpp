
#ifndef HITRUN_HPP
#define HITRUN_HPP

#include "RandomWalk.hpp"

class HitAndRun: public RandomWalk{

    public:
        /**
         * @brief Hit and Run implementation constructor
         * @param r spread hyperparamter
         * @param err error hyperparameter
         * @param thin thin parameter (record every ith value)
         */
        HitAndRun(double r, double err = 1e-6, int thin = 1) : ERR(err), R(r), RandomWalk(thin) {

        }

        /**
         * @brief Generate values from the walk
         * @param num_steps number of steps wanted to take
         * @param x initial starting point
         * @param A polytope matrix
         * @param b polytope matrix
         * @return Matrix
         */
        MatrixXd generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b, int burn) override;

         /**
         * @brief print general type 
         * @return void
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