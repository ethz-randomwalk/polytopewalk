#ifndef LEVSCORE_HPP
#define LEVSCORE_HPP
#include "Common.hpp"

class LeverageScore{
    public:

    LeverageScore (){};
    /**
     * @brief get the Leverage Score approximate calculation
     * @param A polytope matrix (Ax = b)
     * @param W Weight Matrix for slack
     * @param x polytope vector (Ax = b)
     * @param ERR error term
     * @param k last k values have inequality constraint
     * @return Vector
     */
    VectorXd generate(const SparseMatrixXd& A, const SparseMatrixXd& W, const VectorXd& x, const double ERR, const int k);
};

#endif