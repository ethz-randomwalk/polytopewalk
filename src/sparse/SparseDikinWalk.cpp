#include "SparseDikinWalk.hpp"

SparseMatrixXd SparseDikinWalk::generateWeight(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){

    return SparseMatrixXd(VectorXd::Ones(A.cols()).asDiagonal());
}

void SparseDikinWalk::setDistTerm(int d, int n){
    DIST_TERM = (R * R)/d; 
}
