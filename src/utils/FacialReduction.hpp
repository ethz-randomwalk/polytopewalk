#ifndef SPARSE_FR_HPP
#define SPARSE_FR_HPP

#include "SparseLP.hpp"

struct z_res{
    bool found_sol; 
    VectorXd z; 
};

struct fr_res{
    SparseMatrixXd A;
    VectorXd b; 
};

struct res{
    SparseMatrixXd sparse_A;
    VectorXd sparse_b; 
    SparseMatrixXd saved_V; 
    MatrixXd dense_A;
    VectorXd dense_b; 
    MatrixXd Q; 
    VectorXd z1;
};

class FacialReduction {
    public:
        /**
         * @brief Facial Reduction initialization
         * @param max_iter maximum iterations during linear program
         * @param tol tolerance term
         * @param s_max error term
         * @param err_lp error sensitivity for lp calculation
         * @param err_dc error sensitivity for decomposition calculation
         * @return res
         */
        FacialReduction(int max_iter = 3000, double tol = 1e-8, double s_max = 100, double err_lp = 1e-8, double err_dc = 1e-5) : MAX_ITER(max_iter), TOL(tol), S_MAX(s_max), ERR_LP(err_lp),  ERR_DC(err_dc){}
        /**
         * @brief completes facial reduction on Ax = b, x >=_c 0
         * @param A polytope matrix (Ax = b)
         * @param b polytope vector (Ax = b)
         * @param k k values >= 0 constraint
         * @param sparse decision to choose full-dimensional or constraint formulation
         * @return res
         */
        res reduce(SparseMatrixXd A, VectorXd b, int k, bool sparse);
    
    protected:
        /**
         * @brief finds a vector z satisfying A^Ty = [0 z], z in R^n, z >= 0, z != 0, <b, y> = 0
         * @param A polytope matrix (Ax = b)
         * @param b polytope vector (Ax = b)
         * @param k values >= 0 constraint
         * @return z_res
         */
        z_res findZ(const SparseMatrixXd& A, const VectorXd& b, int k);

        /**
         * @brief finds supports with z vector
         * @param z vector
         * @param k values >= 0 constraint
         * @return SparseMatrixXd 
         */
        SparseMatrixXd pickV(const VectorXd& z, int k);

        /**
         * @brief removes redundant constraints in AV
         * @param AV matrix to remove redundant constraints
         * @return SparseMatrixXd 
         */
        SparseMatrixXd pickP(const SparseMatrixXd& AV);

        /**
         * @brief iteratively reduces dimension of the problem using recursion
         * @param A polytope matrix (Ax = b)
         * @param b polytope vector (Ax = b)
         * @param k values >= 0 constraint
         * @return fr_res
         */
        fr_res entireFacialReductionStep(SparseMatrixXd A, VectorXd b, int k);

        /**
         * @brief return global V variable
         */
        SparseMatrixXd savedV;

        /**
         * @brief max iterations on linear program
         */
        const int MAX_ITER;

        /**
         * @brief tolerance parameter
         */
        const double TOL;

        /**
         * @brief tolerance parameter
         */
        const double S_MAX;

         /**
         * @brief LP error parameter
         */
        const double ERR_LP; 

        /**
         * @brief DC error parameter
         */
        const double ERR_DC; 

        /**
         * @brief save last index
         */
        int global_index; 
};

#endif