#include "DenseCenter.hpp"

VectorXd DenseCenter::getInitialPoint(MatrixXd A, VectorXd b){
    glp_prob *lp;
    lp = glp_create_prob();
    glp_term_out(GLP_OFF);
    int amount = 1 + (A.rows() * (A.cols() + 1)); 
    int ia [amount];
    int ja [amount];
    double ar [amount];

    int row_length = A.rows(); 
    int col_length = A.cols() + 1; 

    glp_add_rows(lp, row_length);
    glp_add_cols(lp, col_length);
    glp_set_obj_coef(lp, col_length , 1);
    glp_set_obj_dir(lp, GLP_MAX);

    for(int i = 0; i < b.rows(); i++){
        glp_set_row_bnds(lp, i + 1, GLP_UP, b(i), b(i));
    }
    for(int i = 0; i < col_length; i++){
        glp_set_col_bnds(lp, i + 1, GLP_FR, 0, 0);
    }

    int ind = 1;
    for(int i = 0; i < A.rows(); i++){
        for(int j = 0; j < A.cols(); j++){
            ia[ind] = i + 1;
            ja[ind] = j + 1;
            ar[ind] = A.coeff(i, j); 
            ind ++; 
        }
        ia[ind] = i + 1;
        ja[ind] = A.cols() + 1; 
        ar[ind] = 1.0; 
        ind ++;
    }

    glp_load_matrix(lp, ind-1, ia, ja, ar);
    glp_simplex(lp, NULL);
    double val = glp_get_obj_val(lp); 

    VectorXd ans(A.cols());
    for(int i = 0; i < A.cols(); i++){
        ans.coeffRef(i) = glp_get_col_prim(lp, i + 1);
    }
    glp_delete_prob(lp);
    return ans; 

}