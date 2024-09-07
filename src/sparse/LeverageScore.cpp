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

    SparseMatrixXd L_col = (perm * L0.transpose() * perm);
    SparseMatrix<double,Eigen::RowMajor> L_row = L_col;

    SparseMatrixXd inv(L0.rows(), L0.rows());

    for(int i = 0; i < L_row.outerSize(); i++){
        for(SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(L_row, i); it; ++it){
            int j = it.col();
            double z = (i == j) ? (double)1/D(L_row.outerSize() - 1 - i) : 0;

            for(SparseMatrix<double,Eigen::RowMajor>::InnerIterator it2(L_row, i); it2; ++it2){
                if (it2.col() >= i) break;
                double val = it.col() <= j ? inv.coeffRef(it2.col(), j) : inv.coeffRef(j, it2.col());
                z -= it2.value() * val;
            }
            inv.coeffRef(i, j) = z;
            inv.coeffRef(j, i) = z;
        }
    }
    
    inv = perm * inv.transpose() * perm;

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
