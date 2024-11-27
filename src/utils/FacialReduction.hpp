#ifndef SPARSE_FR_HPP
#define SPARSE_FR_HPP

#include "Common.hpp"
#include "SparseLP.hpp"

/**
* @brief result of Find Z algorithm 
* @param found_sol if the algorithm found a z
* @param z the vector z
*/
struct ZResult{
    bool found_sol; 
    VectorXd z; 
};

/**
* @brief result of Facial Reduction step 
* @param A Ax = b
* @param b Ax = b
* @param savedV PAVv = Pb decomposition
*/
struct FRResult{
    SparseMatrixXd A;
    VectorXd b; 
    SparseMatrixXd savedV;
};

/**
* @brief final output of Facial Reduction algorithm 
* @param sparse_A constrained form Ax = b, x >=_k 0
* @param sparse_b constrained form Ax = b, x >=_k 0
* @param sparse_V PAVv = Pb decomposition
* @param dense_A full-dim form Ax <= b
* @param dense_b full-dim form Ax <= b
* @param Q matrix used to go between forms
* @param z1 vector used to go between forms
*/
struct FROutput{
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
         * @brief initialization for Facial Reduction class
         * @param err_dc error sensitivity for decomposition calculation
         */
        FacialReduction(double err_dc = 1e-5) : ERR_DC(err_dc){}
        /**
         * @brief completes facial reduction on Ax = b, x >=_k 0
         * @param A polytope matrix (Ax = b)
         * @param b polytope vector (Ax = b)
         * @param k k values >= 0 constraint
         * @param sparse decision to choose full-dimensional or constraint formulation
         * @return FROutput
         */
        FROutput reduce(SparseMatrixXd A, VectorXd b, int k, bool sparse);
    
    protected:
        /**
         * @brief finds a vector z satisfying A^Ty = [0 z], z in R^n, z >= 0, z != 0, <b, y> = 0
         * @param A polytope matrix (Ax = b)
         * @param b polytope vector (Ax = b)
         * @param k values >= 0 constraint
         * @return ZResult
         */
        ZResult findZ(const SparseMatrixXd& A, const VectorXd& b, int k);

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
         * @param savedV V in AVv = b
         * @return FRResult
         */
        FRResult entireFacialReductionStep(SparseMatrixXd& A, VectorXd& b, int k, SparseMatrixXd& savedV);

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