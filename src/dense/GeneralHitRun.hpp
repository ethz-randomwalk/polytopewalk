
#ifndef GENHITRUN_HPP
#define GENHITRUN_HPP

#include "RandomWalk.hpp"

class GeneralHitAndRun: public RandomWalk{

    public:
        /**
         * @brief initialization of General Hit and Run class
         * @param r spread hyperparameter
         * @param dist_func distribution
         * @param thin thin parameter (record every ith value)
         * @param err error hyperparameter
         */
        GeneralHitAndRun(double r, function<double(const VectorXd&)> dist_func, int thin = 1,
        double err = 1e-6) : R(r), DIST_FUNC(dist_func), ERR(err), RandomWalk(thin) {}

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
         * @brief distribution type {uniform, normal, log-concave}
         */
        const function<double(const VectorXd&)> DIST_FUNC;

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
         * @param A polytope matrix
         * @param b polytope vector
         * @return double 
         */
        double binarySearch(VectorXd direction, VectorXd& x, const MatrixXd& A, const VectorXd& b);

};

#endif