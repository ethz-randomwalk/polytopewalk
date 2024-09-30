#include "FacialReduction.hpp"


z_res FacialReduction::findZ(const SparseMatrixXd& A, const VectorXd& b, int x_dim){
    z_res ans;
    ans.found_sol = false;
    SparseLP sparse_lp;

    int row_length = A.cols() + 2;
    int col_length = A.rows();
    SparseMatrixXd obj_mat(row_length, col_length);
    VectorXd obj_vec = VectorXd::Zero(col_length);
    obj_vec(obj_vec.rows() - 1) = 1; 

    VectorXd row_bnds = VectorXd::Zero(row_length);
    VectorXd row_rel = VectorXd::Zero(row_length);
    VectorXd col_bnds = VectorXd::Zero(col_length);
    VectorXd col_rel = VectorXd::Zero(col_length); 

    for(int i = 0; i < x_dim; i++){
        row_bnds(i) = 0;
        row_rel(i) = GLP_FX; 
    }
    for(int i = x_dim; i < A.cols(); i++){
        row_bnds(i) = 0;
        row_rel(i) = GLP_LO; 
    }
    row_bnds(A.cols()) = 0;
    row_rel(A.cols()) = GLP_FX;
    row_bnds(A.cols() + 1) = 1;
    row_rel(A.cols() + 1) = GLP_FX; 

    for(int i = 0; i < col_length; i++){
        col_bnds(i) = 0;
        col_rel(i) = GLP_FR; 
    }
    for(int i = 0; i < A.outerSize(); i++){
        for(SparseMatrixXd::InnerIterator it(A, i); it; ++it){
            int row = it.row();
            int col = it.col();
            double val = it.value();

            obj_mat.insert(col, row) = val; 
        }
    }
    for(int i = 0; i < b.rows(); i++){
        obj_mat.insert(A.cols(), i) = b.coeff(i); 
    }
    for(int i = global_index; i < A.cols(); i++){
        for(int j = 0; j < A.rows(); j++){
            double val = A.coeff(j, i); 
            obj_mat.coeffRef(A.cols() + 1, j) = val;
        }

        VectorXd sol = sparse_lp.findOptimalVector(obj_mat, row_bnds, obj_vec, row_rel, col_bnds, col_rel);
        if (sol.cwiseAbs().sum() != 0){
            ans.found_sol = true; 
            ans.z = (A.transpose() * sol);
            return ans;
        }
        global_index++;
    }
    return ans; 
}

SparseMatrixXd FacialReduction::pickV(const VectorXd& z, int x_dim){
    int d = z.rows();
    vector<T> indices;
    for(int i = 0; i < x_dim; i++){
        indices.push_back(T(indices.size(), i, 1)); 
    }
    for(int i = x_dim; i < d; i++){
         if(z(i) < ERR_DC) indices.push_back(T(indices.size(), i, 1)); 
    }
    SparseMatrixXd mat(indices.size(), d);
    mat.setFromTriplets(indices.begin(), indices.end());
    return mat.transpose();
}

SparseMatrixXd FacialReduction::pickP(const SparseMatrixXd& AV){
    SparseQR<SparseMatrixXd, NaturalOrdering<SparseMatrix<double>::StorageIndex>> solver;
    solver.compute(AV.transpose());
    SparseMatrixXd R = solver.matrixR();

    vector<T> indices;
    for (int i = 0; i < min(R.cols(), R.rows()); i++){
        if (abs(R.coeffRef(i, i)) > ERR_DC){
            indices.push_back(T(indices.size(), solver.colsPermutation().indices()(i), 1));
        }
    }
    SparseMatrixXd proj (indices.size(), AV.rows());
    proj.setFromTriplets(indices.begin(), indices.end());
    return proj; 
}

fr_res FacialReduction::entireFacialReductionStep(SparseMatrixXd A, VectorXd b, int x_dim){
    z_res z_ans = findZ(A, b, x_dim);

    if(!z_ans.found_sol){
        fr_res ans;
        ans.A = A;
        ans.b = b; 
        return ans; 
    }
    SparseMatrixXd V = pickV(z_ans.z, x_dim);
    savedV = savedV * V; 
    SparseMatrixXd AV = A * V;
    SparseMatrixXd P = pickP(AV);
    A = P * AV;
    b = P * b; 
    return entireFacialReductionStep(A, b, x_dim);
}

res FacialReduction::reduce(SparseMatrixXd A, VectorXd b, int k, bool sparse){
    int x_dim = A.cols() - k; 
    savedV = SparseMatrixXd(VectorXd::Ones(A.cols()).asDiagonal());
    global_index = x_dim; 
    //remove dependent rows
    SparseMatrixXd P = pickP(A);
    A = P * A; 
    b = P * b; 
    fr_res result = entireFacialReductionStep(A, b, x_dim);
    res final_res; 

    final_res.sparse_A = result.A;
    final_res.sparse_b = result.b;
    final_res.saved_V = savedV;

    if(!sparse){
        HouseholderQR <MatrixXd> qr(result.A.cols(), result.A.rows());
        qr.compute(MatrixXd(result.A.transpose()));
        MatrixXd Q = qr.householderQ();
        MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
        int d = R.rows();
        int n = R.cols();

        MatrixXd newR = R.block(0, 0, R.cols(), R.cols());
        VectorXd z1 = newR.transpose().inverse() * result.b;

        MatrixXd Q1 = Q.block(0, 0, Q.rows(), n);
        MatrixXd Q2 = Q.block(0, n, Q.rows(), d - n);
        MatrixXd reduced_A = -1 * Q2.block(x_dim, 0, d - x_dim, d - n);
        VectorXd reduced_b = (Q1 * z1).tail(d - x_dim);

        final_res.dense_A = reduced_A;
        final_res.dense_b = reduced_b;
        
        final_res.z1 = z1;
        final_res.Q = Q;
    }
    return final_res;
}

