
#ifndef DIKINWALK_HPP
#define DIKINWALK_HPP

#include "BarrierWalk.hpp"

class DikinWalk: public BarrierWalk{

    public:
        /**
         * @brief constructor for Dikin Walk class
         * @param r spread parameter
         * @param thin thin parameter
         */
        DikinWalk(double r, int thin = 1) : BarrierWalk(r, thin){}

        /**
         * @brief print dikin
         * @return void
         */
        void printType() override;

        /**
         * @brief returns weight for DikinWalk (Identity Matrix)
         * @param x point
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @return void
         */
        void generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd&b) override;

    protected:

        /**
         * @brief set Dist Term for Dikin Walk
         * @param d (dimension)
         * @param n (number of constraints)
         * @return void
         */
        void setDistTerm(int d, int n) override;

};


#endif