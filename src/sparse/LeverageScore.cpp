#include "LeverageScore.hpp"

VectorXd LeverageScore::generate(const SparseMatrixXd& A, const SparseMatrixXd& W, const VectorXd& x, const double ERR, const int k){
    
    SparseMatrixXd S_inv(x.rows(), x.rows());
    for(int i = x.rows() - k; i < x.rows(); i++){
        S_inv.coeffRef(i, i) = 1/x(i);
    }
    SparseMatrixXd G_sqrt = S_inv * W;
    for(int i = 0; i < x.rows() - k; i++){
        G_sqrt.coeffRef(i, i) = ERR;
    }

    SparseMatrixXd G_inv_sqrt = SparseMatrixXd(VectorXd(G_sqrt.diagonal()).cwiseInverse().asDiagonal());

    SparseMatrixXd hess = A * (G_inv_sqrt * G_inv_sqrt)* A.transpose();
    SimplicialLDLT<SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> cholesky;
    cholesky.analyzePattern(hess);
    cholesky.factorize(hess);
    SparseMatrixXd L0 = cholesky.matrixL();
    VectorXd D = cholesky.vectorD();
    SparseMatrixXd D_sqrt (L0.rows(), L0.cols());

    for(int i = 0; i < L0.rows(); i++){
        D_sqrt.coeffRef(i, i) = sqrt(D(i));
    }
    SparseMatrixXd L = L0 * D_sqrt;
    map<int, vector<int>> sparsity_ref; 
    // row : list of col indices
    for(int p = 0; p < L.outerSize(); p++){
        for(SparseMatrixXd::InnerIterator it(L, p); it; ++it){
            sparsity_ref[it.col()].push_back(it.row());
        }
    }

    SparseMatrixXd inv(L.rows(), L.rows());
    map<int, vector<pair<int, int>>> row_major;
    map<int, vector<pair<int, int>>> col_major;

    for(int i = L.rows() - 1; i >= 0; i--){
        for(int v = sparsity_ref[i].size() - 1; v >= 0; v--){
            int j = sparsity_ref[i][v];
            double z = (i == j) ? (double)1/D(i) : 0;

            for(auto pair : col_major[j]){
                int first = pair.first;
                int second = pair.second;
                double term = inv.coeff(first, second);
                int p = min(first, second);

                z -= term * L0.coeff(p, i);
            }

            for(auto pair : row_major[j]){
                int first = pair.first;
                int second = pair.second;
                double term = inv.coeff(first, second);
                int p = max(first, second);
                z -= term * L0.coeff(p, i);
            }

            col_major[j].push_back(make_pair(i, j));
            if (i != j){
                row_major[i].push_back(make_pair(i, j));
            }

            inv.coeffRef(i, j) = z;
            inv.coeffRef(j, i) = z; 
        }
    }

    SparseMatrixXd AG_inv_sqrt = A * G_inv_sqrt;
    SparseMatrixXd P = inv * AG_inv_sqrt;
    SparseMatrixXd result (AG_inv_sqrt.cols(), AG_inv_sqrt.cols());
    for(int i = 0; i < AG_inv_sqrt.cols(); i++){
        double val = AG_inv_sqrt.col(i).dot(P.col(i));
        result.coeffRef(i, i) = val;
    }

    SparseMatrixXd I = SparseMatrixXd(x.rows(), x.rows());
    for(int i = 0; i < x.rows() - k; i++){
        I.coeffRef(i, i) = result.coeffRef(i, i);
    }
    for(int i = x.rows() - k; i < x.rows(); i++){
        I.coeffRef(i, i) = 1; 
    }
    result = I - result;

    return result.diagonal();
}
