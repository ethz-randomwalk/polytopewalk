#include "SparseCenter.hpp"

VectorXd SparseCenter::getInitialPoint(SparseMatrixXd& A, VectorXd& b, int k){

    // Solve the linear program
    // max delta 
    // s.t. A x = b
    // and x>= delta, on the last k coordinates 
    SparseLP sparse_lp;
    int row_length = A.rows() + k;
    int col_length = A.cols() + 1; 

    SparseMatrixXd obj_mat(row_length, col_length);
    VectorXd obj_vec = VectorXd::Zero(col_length);
    obj_vec(obj_vec.rows() - 1) = 1; 

    VectorXd row_bnds = VectorXd::Zero(row_length);
    VectorXd row_rel = VectorXd::Zero(row_length);
    VectorXd col_bnds = VectorXd::Zero(col_length);
    VectorXd col_rel = VectorXd::Zero(col_length); 
    vector<T> coefficients;  

    for(int i = 0; i < b.rows(); i++){
        row_bnds(i) = b(i);
        row_rel(i) = GLP_FX; 
    }
    for(int i = b.rows(); i < row_length; i++){
        row_rel(i) = GLP_LO; 
    }
    for(int i = 0; i < col_length - k - 1; i++){
        col_rel(i) = GLP_FR; 
    }
    for(int i = col_length - k - 1; i < col_length; i++){
        col_rel(i) = GLP_LO; 
    }

    for(int i = 0; i < A.outerSize(); i++){
        for(SparseMatrixXd::InnerIterator it(A, i); it; ++it){
            int row = it.row();
            int col = it.col();
            double val = it.value();
            coefficients.push_back(T(row, col, val));
        }
    }
    for(int i = 0; i < k; i++){
        int row_val = A.rows() + i; 
        int col_val = A.cols() - k + i; 
        coefficients.push_back(T(row_val, col_val, 1));
        coefficients.push_back(T(row_val, A.cols(), -1));
    }
    obj_mat.setFromTriplets(coefficients.begin(), coefficients.end());

    // call the lp solver
    VectorXd sol = sparse_lp.findOptimalVector(obj_mat, row_bnds, obj_vec, row_rel, col_bnds, col_rel); 

    // retrieve x 
    VectorXd ans = VectorXd::Zero(sol.rows() - 1);
    for(int i = 0; i < ans.rows(); i++){
        ans(i) = sol(i);
    }
    return ans;
    
}