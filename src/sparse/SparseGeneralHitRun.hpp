#ifndef CONSGENHITRUN_HPP
#define CONSGENHITRUN_HPP

#include "SparseRandomWalk.hpp"

class SparseGeneralHitAndRun : public SparseRandomWalk{
    public:
        /**
         * @brief initialization of Sparse Hit and Run class
         * @param r spread parameter
         * @param dist_func distribution function
         * @param thin thin parameter
         * @param err error constant
         */
        SparseGeneralHitAndRun(double r, function<double(const VectorXd&)> dist_func, 
                int thin = 1, double err = 1e-6) : R(r), DIST_FUNC(dist_func), SparseRandomWalk(thin, err) {}

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
         * @brief distribution type {uniform, normal, log-concave}
         */
        const function<double(const VectorXd&)> DIST_FUNC;


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