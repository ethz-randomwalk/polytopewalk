#include "SparseJohnWalk.hpp"
#include "LeverageScore.hpp"

SparseMatrixXd SparseJohnWalk::generateWeight(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    double d = A.cols() - A.rows();
    double n = k; 
    double alpha = 1 - 1/(log2(2.0 * n / d));
    double beta = (double)d / (2.0 * n);

    if (w_i.coeffRef(0) == -1 || w_i.rows() != x.rows()){
        w_i = VectorXd::Ones(x.rows());
    }

    LeverageScore L;
    VectorXd beta_ones = beta * VectorXd::Ones(x.rows());

    for(int i = 0; i < w_i.rows() - k; i++){
        w_i(i) = 0;
        beta_ones.coeffRef(i) = 0;
    }
    VectorXd next_weight = w_i;

    // fixed point iteration
    for(int i = 0; i < MAX_ITER; i++){
        w_i = next_weight;
        SparseMatrixXd W (w_i.rows(), w_i.rows());
        for(int j = x.rows() - k; j < x.rows(); j++){
            W.coeffRef(j, j) = pow(w_i(j), alpha * 0.5);
        }
        VectorXd score  = L.generate(A, W, x, ERR, k);
        next_weight = 0.5 * (w_i + score + beta_ones).cwiseMax(beta_ones);

        if ((w_i - next_weight).cwiseAbs().maxCoeff() < LIM){
            break;
        }
    }

    return SparseMatrixXd(w_i.asDiagonal());

}

void SparseJohnWalk::setDistTerm(int d, int n){
    w_i = VectorXd::Ones(d);
    DIST_TERM = (R * R)/pow(d, 1.5);
}