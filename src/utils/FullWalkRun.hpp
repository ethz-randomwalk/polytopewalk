#include "FacialReduction.hpp"
#include "DenseCenter.hpp"
#include "SparseCenter.hpp"
#include "dense/DikinWalk.hpp"
#include "dense/DikinLSWalk.hpp"
#include "dense/JohnWalk.hpp"
#include "dense/VaidyaWalk.hpp"
#include "dense/HitRun.hpp"
#include "dense/BallWalk.hpp"
#include "sparse/SparseDikinWalk.hpp"
#include "sparse/SparseDikinLSWalk.hpp"
#include "sparse/SparseJohnWalk.hpp"
#include "sparse/SparseVaidyaWalk.hpp"
#include "sparse/SparseBallWalk.hpp"
#include "sparse/SparseHitRun.hpp"


/**
 * @brief runs full preprocessing, walk, and post-processing steps in dense formulation
 * @param niter number of steps
 * @param A polytope matrix (Ax = b)
 * @param b polytope vector (Ax = b)
 * @param k values >= 0 constraint
 * @param walk dense random walk implementation
 * @param fr facial reduction algorithm
 * @param init initialization algorithm 
 * @param burnin how many to exclude
 * @param thin thinning parameter
 * @param seed seed for reproducibility
 * @return (niter - burnin)//thin by d (dimension of x) matrix
 */
MatrixXd denseFullWalkRun(int niter, SparseMatrixXd A, VectorXd b, int k, RandomWalk* walk, FacialReduction* fr, DenseCenter* init, int burnin = 0, int thin = 1, int seed = -1){
    if (k <= 0 || k > A.cols()) {
        throw std::invalid_argument("Parameter k must be between 1 and the number of columns in A.");
    }
    
    FROutput fr_result = fr->reduce(A, b, k, false);
    VectorXd x = init->getInitialPoint(fr_result.dense_A, fr_result.dense_b);
    MatrixXd steps = walk->generateCompleteWalk(niter, x, fr_result.dense_A, fr_result.dense_b, burnin, thin, seed);
    MatrixXd res(steps.rows(), A.cols());
    for(int i = 0; i < steps.rows(); i++){
        VectorXd val (steps.cols() + fr_result.z1.rows());
        VectorXd row = steps.row(i);
        val << fr_result.z1, row;
        res.row(i) = (fr_result.Q * val).head(A.cols());
    }
    return res; 
}

/**
 * @brief runs full preprocessing, walk, and post-processing steps in sparse formulation
 * @param niter number of steps
 * @param A polytope matrix (Ax <= b)
 * @param b polytope vector (Ax <= b)
 * @param k last k coordinates >= 0
 * @param walk sparse random walk implementation
 * @param fr facial reduction algorithm
 * @param init initialization algorithm 
 * @param burnin how many to exclude
 * @param thin thinning parameter
 * @param seed seed for reproducibility
 * @return (niter - burnin)//thin by d (dimension of x) matrix
 */
MatrixXd sparseFullWalkRun(int niter, SparseMatrixXd A, VectorXd b, int k, SparseRandomWalk* walk, FacialReduction* fr, SparseCenter* init, int burnin = 0, int thin = 1, int seed = -1){
    if (k <= 0 || k > A.cols()) {
        throw std::invalid_argument("Parameter k must be between 1 and the number of columns in A.");
    }
    
    FROutput fr_result = fr->reduce(A, b, k, true);
    int new_k = fr_result.sparse_A.rows() - (A.rows() - k);
    VectorXd x = init->getInitialPoint(fr_result.sparse_A, fr_result.sparse_b, new_k);
    MatrixXd steps = walk->generateCompleteWalk(niter, x, fr_result.sparse_A, fr_result.sparse_b, new_k, burnin, thin, seed);
    MatrixXd res(steps.rows(), A.cols());
    for(int i = 0; i < steps.rows(); i++){
        res.row(i) = fr_result.saved_V * steps.row(i).transpose();
    }
    return res; 
}