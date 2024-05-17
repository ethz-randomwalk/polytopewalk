#include "SparseVaidyaWalk.hpp"
#include "LeverageScore.hpp"

SparseMatrixXd SparseVaidyaWalk::generateWeight(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    LeverageScore L;
    SparseMatrixXd W (x.rows(), x.rows());
    for(int i = x.rows() - k; i < x.rows(); i++){
        W.coeffRef(i, i) = 1;
    }
    VectorXd weights = L.generate(A, W, x, ERR, k);
    for (int i = weights.rows() - k; i < weights.rows(); i++){
        weights(i) += ((double)(A.cols() - A.rows())/k);
    }
    return SparseMatrixXd(weights.asDiagonal());
}

void SparseVaidyaWalk::setDistTerm(int d, int n){
    DIST_TERM = (R * R)/sqrt(n * d);
}