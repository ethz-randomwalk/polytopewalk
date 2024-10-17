#ifndef CPF_HPP
#define CPF_HPP

#include "Common.hpp"


class DenseCenter {
    public:

        /**
         * @brief Center Algorithm
         * @return res
         */
        DenseCenter(){};

        /**
         * @brief Finds analytical center Ax <= b
         * @param A polytope matrix (Ax <= b)
         * @param b polytope vector (Ax <= b)
         * @return VectorXd 
         */
        VectorXd getInitialPoint(MatrixXd& A, VectorXd& b);

};

#endif