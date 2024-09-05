#include "LeverageScore.hpp"
#include <chrono>

VectorXd LeverageScore::generate(const SparseMatrixXd& A, const SparseMatrixXd& W, const VectorXd& x, const double ERR, const int k){
    

    auto overall_start = chrono::high_resolution_clock::now();
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
    auto start = chrono::high_resolution_clock::now();

    SimplicialLDLT<SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> cholesky;
    cholesky.analyzePattern(hess);
    cholesky.factorize(hess);
    SparseMatrixXd L0 = cholesky.matrixL();
    VectorXd D = cholesky.vectorD();
    SparseMatrixXd D_sqrt (L0.rows(), L0.cols());
    auto stop = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "CHOLESKY DURATION" << endl;
    cout << duration.count() << endl;

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
    start = chrono::high_resolution_clock::now();
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
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
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

    auto overall_stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(overall_stop - overall_start);
    cout << "OVERALL DURATION" << endl;
    cout << duration.count() << endl;
}
