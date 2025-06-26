
#ifndef BALLWALK_HPP
#define BALLWALK_HPP

#include "RandomWalk.hpp"

class BallWalk: public RandomWalk{
    

    public:

        /**
         * @brief initialization of Ball Walk class
         * @param r spread parameter
         */
        BallWalk(double r) : R(r){}

        /**
         * @brief generate values from Ball Walk
         * @param niter number of steps wanted to take
         * @param init initial starting point
         * @param A polytope matrixd (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @param burnin number of initial steps to cut
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
         * @brief spread parameter
         */
        const double R;


};

#endif