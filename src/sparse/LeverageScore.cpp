#include "LeverageScore.hpp"

VectorXd LeverageScore::generate(const SparseMatrixXd& A, const SparseMatrixXd& W, const VectorXd& x, const double ERR, const int k){
    
    VectorXd S_inv(x.rows());
    for(int i = x.rows() - k; i < x.rows(); i++){
        S_inv.coeffRef(i) = 1/x(i);
    }
    VectorXd G_sqrt = W * S_inv;
    for(int i = 0; i < x.rows() - k; i++){
        G_sqrt.coeffRef(i) = ERR;
    }

    SparseMatrixXd G_inv_sqrt = SparseMatrixXd(G_sqrt.cwiseInverse().asDiagonal());
    SparseMatrixXd AG_inv_sqrt = A * G_inv_sqrt;

    SparseMatrixXd hess = AG_inv_sqrt * AG_inv_sqrt.transpose();

    SimplicialLDLT<SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> cholesky;
    cholesky.analyzePattern(hess);
    cholesky.factorize(hess);
    SparseMatrixXd L0 = cholesky.matrixL();
    VectorXd D = cholesky.vectorD();

    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm (L0.rows());
    for(int i = 0; i < L0.rows(); i++){
        perm.indices()(i) = perm.rows() - 1 - i;
    }

    SparseMatrixXd L_col = perm * L0 * perm;

    SparseMatrixXd inv(L0.rows(), L0.rows());

    VectorXd nnz (L0.rows());
    for(int i = 0; i < L_col.rows(); i++){
        nnz(i) = L_col.col(i).nonZeros();
    }
    inv.reserve(nnz);

    for(int i = 0; i < L_col.outerSize(); i++){
        for(SparseMatrixXd::InnerIterator it(L_col, i); it; ++it){
            int j = it.row();
            double z = (i == j) ? (double)1/D(L_col.outerSize() - 1 - i) : 0;

            for(SparseMatrixXd::InnerIterator it2(L_col, i); it2; ++it2){
                if (it2.row() >= i) break;
                double val = it2.row() <= j ? inv.coeffRef(it2.row(), j) : inv.coeffRef(j, it2.row());
                z -= it2.value() * val;
            }
            if (i <= j) inv.insert(i, j) = z;
            else inv.insert(j, i) = z;
        }
    }
    inv = perm * inv * perm;

    SparseMatrixXd P = inv.selfadjointView<Lower>() * AG_inv_sqrt;
    VectorXd result (AG_inv_sqrt.cols());
    for(int i = 0; i < AG_inv_sqrt.cols(); i++){
        double val = AG_inv_sqrt.col(i).dot(P.col(i));
        result.coeffRef(i) = val;
    }

    for(int i = 0; i < x.rows() - k; i++){
        result(i) = 0;
    }

    for(int i = x.rows() - k; i < x.rows(); i++){
        result(i) = 1 - result(i);
    }

    return result; 
}