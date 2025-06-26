
#ifndef VAIDYAWALK_HPP
#define VAIDYAWALK_HPP

#include "BarrierWalk.hpp"

class VaidyaWalk: public BarrierWalk{

    public:
        /**
         * @brief initialization of Sparse Vaidya Walk class
         * @param r spread parameter
         */
        VaidyaWalk(double r) : BarrierWalk(r){}  

        /**
         * @brief print general type 
         */
        void printType() override;

        /**
         * @brief returns weight for Vaidya Walk (leverage score calculation)
         * @param x center vector
         * @param A polytope matrix
         * @param b polytope vector
         * @returns Vectorr
         */
        VectorXd generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b) override;
    
    protected:

        /**
         * @brief set distribution constant
         * @param d (dimension)
         * @param n (number of constraints)
         */
        void setDistTerm(int d, int n) override;

};

#endif