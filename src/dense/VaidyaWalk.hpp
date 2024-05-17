
#ifndef VAIDYAWALK_HPP
#define VAIDYAWALK_HPP

#include "BarrierWalk.hpp"

class VaidyaWalk: public BarrierWalk{

    public:
        /**
         * @brief constructor for Vaidya Walk class
         * @param r spread parameter
         * @param thin thin constant
         */
        VaidyaWalk(double r, int thin = 1) : BarrierWalk(r, thin){}  

        /**
         * @brief print general type 
         * @return void
         */
        void printType() override;

        /**
         * @brief returns weight for Vaidya Walk (leverage score calculation)
         * @param x center vector
         * @param A polytope matrix
         * @param b polytope vector
         * @return void
         */
        void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b) override;
    
    protected:

        /**
         * @brief set Dist Term for Vaidya Walk
         * @param d (dimension)
         * @param n (number of constraints)
         * @return void
         */
        void setDistTerm(int d, int n) override;

};

#endif