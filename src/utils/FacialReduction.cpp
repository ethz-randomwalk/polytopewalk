#include "FacialReduction.hpp"


z_res FacialReduction::findZ(const SparseMatrixXd& A, const VectorXd& b, int x_dim){
    z_res ans;
    ans.found_sol = false;

    glp_prob *lp;
    glp_term_out(GLP_OFF);
    lp = glp_create_prob();
    int amount = 1 + A.nonZeros() + 2 * b.rows();
    int ia [amount];
    int ja [amount];
    double ar [amount];

    int row_length = A.cols() + 2;
    int col_length = A.rows();

    glp_add_rows(lp, row_length);
    glp_add_cols(lp, col_length);
    glp_set_obj_coef(lp, col_length, 1);

    for(int i = 1; i <= col_length; i++){
        glp_set_col_bnds(lp, i, GLP_FR, 0, 0);
    }

    for(int i = 0; i < x_dim; i++){
        glp_set_row_bnds(lp, i + 1, GLP_FX, 0, 0);
    }
    
    for(int i = x_dim; i < A.cols(); i++){
        glp_set_row_bnds(lp, i + 1, GLP_LO, 0, 0);
    }
    glp_set_row_bnds(lp, A.cols() + 1, GLP_FX, 0, 0);
    glp_set_row_bnds(lp, A.cols() + 2, GLP_FX, 1, 1);

    int ind = 1;
    for(int i = 0; i < A.outerSize(); i++){
        for(SparseMatrixXd::InnerIterator it(A, i); it; ++it){
            int row = it.row();
            int col = it.col();
            double val = it.value();

            ia[ind] = col + 1;
            ja[ind] = row + 1; 
            ar[ind] = val; 
            ind ++; 
        }
    }
    for(int i = 0; i < b.rows(); i++){
        ia[ind] = A.cols() + 1;
        ja[ind] = i + 1;
        ar[ind] = b.coeff(i); 
        ind++;
    }

    int saved_ind = ind;
    for(int i = global_index; i < A.outerSize(); i++){
        ind = saved_ind;
        for(SparseMatrixXd::InnerIterator it(A, i); it; ++it){
            int row = it.row();
            int col = it.col();
            double val = it.value();
            ia[ind] = A.cols() + 2; 
            ja[ind] = row + 1; 
            ar[ind] = val; 
            ind ++; 
        }

        glp_load_matrix(lp, ind - 1, ia, ja, ar);
        glp_simplex(lp, NULL);
        VectorXd sol(col_length);
        for(int i = 0; i < col_length; i++){
            sol.coeffRef(i) = glp_get_col_prim(lp, i + 1);
        }
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

