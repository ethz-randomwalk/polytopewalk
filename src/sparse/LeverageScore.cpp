#include "LeverageScore.hpp"
#include <chrono>

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

    auto start = chrono::high_resolution_clock::now();
    SparseMatrixXd perm (L0.rows(), L0.rows());
    for(int i = 0; i < perm.rows(); i++){
        perm.coeffRef(i, perm.rows() - 1 - i) = 1; 
    }

    SparseMatrixXd L_col = (perm * L0 * perm);

    SparseMatrixXd inv(L0.rows(), L0.rows());

    for(int i = 0; i < L_col.outerSize(); i++){
        for(SparseMatrixXd::InnerIterator it(L_col, i); it; ++it){
            int j = it.row();
            double z = (i == j) ? (double)1/D(L_col.outerSize() - 1 - i) : 0;

            for(SparseMatrixXd::InnerIterator it2(L_col, i); it2; ++it2){
                if (it2.row() >= i) break;
                double val = it.row() <= j ? inv.coeffRef(it2.row(), j) : inv.coeffRef(j, it2.row());
                z -= it2.value() * val;
            }
            inv.coeffRef(i, j) = z;
            inv.coeffRef(j, i) = z;
        }
    }
    
    inv = perm * inv * perm;

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "FOR LOOP DURATION" << endl;
    cout << duration.count() << endl;


    SparseMatrixXd P = inv * AG_inv_sqrt;
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