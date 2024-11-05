#include "SparseDikinLSWalk.hpp"
SparseMatrixXd SparseDikinLSWalk::generateWeight(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    LeverageScore L;

    double d = A.cols() - A.rows();
    double n = k; 
    double q = 2.0 * (1.0 + log(n));
    double alpha = 1.0 - 2.0/q; 

    if (w_i.coeffRef(0) == -1 || w_i.rows() != x.rows()){
        w_i = VectorXd::Ones(x.rows());
    }

    // term1 is all 1 vect on the first k coordinates
    VectorXd term1 = (alpha) * VectorXd::Ones(x.rows());
    VectorXd errors = ERR * VectorXd::Ones(x.rows());

    for(int i = 0; i < w_i.rows() - k; i++){
        w_i(i) = 0;
        term1(i) = 0;
    }
    // gradient descent to optimize the LS barrier
    for(int i = 0; i < MAX_ITER; i++){
        SparseMatrixXd W (x.rows(), x.rows());
        VectorXd term2a = VectorXd::Zero(x.rows());
        for(int j = x.rows() - k; j < x.rows(); j++){
            // term2a = alpha / w
            term2a(j) = (double)alpha/w_i(j);
            W.coeffRef(j, j) = pow(w_i(j), alpha * 0.5);
        }

        // term2b is leverage score
        VectorXd term2b = L.generate(A, W, x, ERR, k);
        // term2 is gradient log det
        // which is the ratio between leverage score and w
        VectorXd term2 = term2a.cwiseProduct(term2b); 
        VectorXd grad = term1 - term2;
        if (grad.norm() < G_LIM){
            break; 
        }
        w_i = (w_i - STEP_SIZE * grad);
        for(int j = x.rows() - k; j < x.rows(); j++){
            w_i(j) = max(w_i(j), ERR);
        }
    }
    return SparseMatrixXd(w_i.asDiagonal());

}

void SparseDikinLSWalk::setDistTerm(int d, int n){
    w_i = VectorXd::Ones(d);
    double q = 2.0 * (1.0 + log(n));
    double term = (1.0 + q) * (1.0 + q * q);
    DIST_TERM = (R * R)/term;
}