#include "DenseCenter.hpp"

VectorXd DenseCenter::getInitialPoint(MatrixXd& A, VectorXd& b){
    // Solve the linear program
    // max delta 
    // s.t. A x + delta * 1 <= b
    glp_prob *lp;
    lp = glp_create_prob();
    glp_term_out(GLP_OFF);
    int amount = 1 + (A.rows() * (A.cols() + 1)); 
    vector<int> ia (amount);
    vector<int> ja (amount);
    vector <double> ar (amount);

    int row_length = A.rows(); 
    // delta is stored at the last column
    int col_length = A.cols() + 1; 

    glp_add_rows(lp, row_length);
    glp_add_cols(lp, col_length);
    // maximize delta * 1
    glp_set_obj_coef(lp, col_length , 1);
    glp_set_obj_dir(lp, GLP_MAX);

    for(int i = 0; i < b.rows(); i++){
        glp_set_row_bnds(lp, i + 1, GLP_UP, b(i), b(i));
    }
    for(int i = 0; i < col_length - 1; i++){
        glp_set_col_bnds(lp, i + 1, GLP_FR, 0, 0);
    }
    glp_set_col_bnds(lp, col_length, GLP_LO, 0, 0);

    int ind = 1;
    for(int i = 0; i < A.rows(); i++){
        // for A * x
        for(int j = 0; j < A.cols(); j++){
            ia[ind] = i + 1;
            ja[ind] = j + 1;
            ar[ind] = A.coeff(i, j); 
            ind ++; 
        }
        // for + delta * 1
        ia[ind] = i + 1;
        ja[ind] = A.cols() + 1; 
        ar[ind] = 1.0; 
        ind ++;
    }

    glp_load_matrix(lp, ind-1, ia.data(), ja.data(), ar.data());
    glp_simplex(lp, NULL);
    double val = glp_get_obj_val(lp); 

    // retrieve x
    VectorXd ans(A.cols());
    for(int i = 0; i < A.cols(); i++){
        ans.coeffRef(i) = glp_get_col_prim(lp, i + 1);
    }
    glp_delete_prob(lp);
    return ans; 

}