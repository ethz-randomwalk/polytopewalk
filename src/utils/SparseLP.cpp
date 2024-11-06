#include "SparseLP.hpp"

VectorXd SparseLP::findOptimalVector(SparseMatrixXd& A, VectorXd& b, VectorXd& c, VectorXd& row_rel, VectorXd& col_cons, VectorXd& col_rel){

    glp_prob *lp;
    glp_term_out(GLP_OFF);
    lp = glp_create_prob();
    int amount = 1 + A.nonZeros();

    vector<int> ia (amount);
    vector<int> ja (amount);
    vector<double> ar (amount);

    int row_length = A.rows();
    int col_length = A.cols(); 
    glp_add_rows(lp, row_length);
    glp_add_cols(lp, col_length);
    for(int i = 0; i < col_length; i++){
        glp_set_obj_coef(lp, i + 1, c(i));
    }
    glp_set_obj_dir(lp, GLP_MAX);
    for(int i = 0; i < row_length; i++){
        glp_set_row_bnds(lp, i + 1, row_rel(i), b(i), b(i));
    }
    for(int i = 0; i < col_length; i++){
        glp_set_col_bnds(lp, i + 1, col_rel(i), col_cons(i), col_cons(i));
    }
    int ind = 1; 
    for(int i = 0; i < A.outerSize(); i++){
        for(SparseMatrixXd::InnerIterator it(A, i); it; ++it){
            int row = it.row();
            int col = it.col();
            double val = it.value();

            ia[ind] = row + 1;
            ja[ind] = col + 1; 
            ar[ind] = val; 

            ind ++; 
        }
    }

    glp_load_matrix(lp, amount-1, ia.data(), ja.data(), ar.data());
    glp_simplex(lp, NULL);
    VectorXd ans(A.cols());
    for(int i = 0; i < A.cols(); i++){
        ans.coeffRef(i) = glp_get_col_prim(lp, i + 1);
    }
    glp_delete_prob(lp);
    return ans;


}