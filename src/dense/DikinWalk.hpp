
#ifndef DIKINWALK_HPP
#define DIKINWALK_HPP

#include "BarrierWalk.hpp"

class DikinWalk: public BarrierWalk{

    public:
        /**
         * @brief initialization of Dikin Walk class
         * @param r spread parameter
         */
        DikinWalk(double r) : BarrierWalk(r){}

        /**
         * @brief print dikin
         */
        void printType() override;

        /**
         * @brief returns weight for DikinWalk (Identity Matrix)
         * @param x point
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @return Vector
         */
        VectorXd generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd&b) override;

    protected:

        /**
         * @brief set distribution constant
         * @param d (dimension)
         * @param n (number of constraints)
         */
        void setDistTerm(int d, int n) override;

};


#endif