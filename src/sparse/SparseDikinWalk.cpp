#include "SparseDikinWalk.hpp"

VectorXd SparseDikinWalk::generateWeight(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){

    return VectorXd::Ones(A.cols());
}

void SparseDikinWalk::setDistTerm(int d, int n){
    DIST_TERM = (R * R)/d; 
}
