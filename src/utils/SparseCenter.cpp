#include "SparseCenter.hpp"

VectorXd SparseCenter::getInitialPoint(SparseMatrixXd A, VectorXd b, int k){
    glp_prob *lp;
    glp_term_out(GLP_OFF);
    lp = glp_create_prob();
    int amount = 1 + A.nonZeros() + 2 * k;
    int ia [amount];
    int ja [amount];
    double ar [amount];

    int row_length = A.rows() + k;
    int col_length = A.cols() + 1; 

    glp_add_rows(lp, row_length);
    glp_add_cols(lp, col_length);
    glp_set_obj_coef(lp, col_length, 1);
    glp_set_obj_dir(lp, GLP_MAX);

    for(int i = 0; i < b.rows(); i++){
        glp_set_row_bnds(lp, i + 1, GLP_FX, b(i), b(i));
    }
    for(int i = b.rows(); i < row_length; i++){
        glp_set_row_bnds(lp, i + 1, GLP_LO, 0, 0);
    }

    for(int i = 0; i < col_length - k - 1; i++){
        glp_set_col_bnds(lp, i + 1, GLP_FR, 0, 0);
    }
    for(int i = col_length - k - 1; i < col_length; i++){
        glp_set_col_bnds(lp, i + 1, GLP_LO, 0, 0);
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
    for(int i = 0; i < k; i++){
        int row_val = A.rows() + 1 + i; 
        int col_val = A.cols() - k + i + 1; 
        ia[ind] = row_val;
        ja[ind] = col_val; 
        ar[ind] = 1.0; 
        ia[ind+1] = row_val;  
        ja[ind+1] = A.cols() + 1; 
        ar[ind+1] = -1.0;
        ind += 2;
    }

    glp_load_matrix(lp, amount-1, ia, ja, ar);
    glp_simplex(lp, NULL);
    double val = glp_get_obj_val(lp); 

    VectorXd ans(A.cols());
    for(int i = 0; i < A.cols(); i++){
        ans.coeffRef(i) = glp_get_col_prim(lp, i + 1);
    }
    glp_delete_prob(lp);
    return ans;
    
}