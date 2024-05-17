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
 * @param A polytope matrix (Ax = b)
 * @param b polytope vector (Ax = b)
 * @param k values >= 0 constraint
 * @param num_sim number of steps
 * @param walk dense random walk implementation
 * @param init initialization algorithm 
 * @return Matrix
 */
MatrixXd denseFullWalkRun(SparseMatrixXd A, VectorXd b, int k, int num_sim, RandomWalk* walk, FacialReduction* fr, DenseCenter* init, int burn = 0){
    res fr_result = fr->reduce(A, b, k, false);
    VectorXd x = init->getInitialPoint(fr_result.dense_A, fr_result.dense_b);
    MatrixXd steps = walk->generateCompleteWalk(num_sim, x, fr_result.dense_A, fr_result.dense_b, burn);
    MatrixXd res(num_sim, A.cols());
    for(int i = 0; i < num_sim; i++){
        VectorXd val (steps.cols() + fr_result.z1.rows());
        VectorXd row = steps.row(i);
        val << fr_result.z1, row;
        res.row(i) = (fr_result.Q * val).head(A.cols());
    }
    return res; 
}

/**
 * @brief runs full preprocessing, walk, and post-processing steps in sparse formulation
 * @param A polytope matrix (Ax <= b)
 * @param b polytope vector (Ax =< b)
 * @param num_sim number of steps
 * @param walk sparse random walk implementation
 * @param init initialization algorithm 
 * @return Matrix
 */
MatrixXd sparseFullWalkRun(SparseMatrixXd A, VectorXd b, int k, int num_sim, SparseRandomWalk* walk, FacialReduction* fr, SparseCenter* init, int burn = 0){
    res fr_result = fr->reduce(A, b, k, true);
    int new_k = fr_result.sparse_A.rows() - (A.rows() - k);
    VectorXd x = init->getInitialPoint(fr_result.sparse_A, fr_result.sparse_b, new_k);
    MatrixXd steps = walk->generateCompleteWalk(num_sim, x, fr_result.sparse_A, fr_result.sparse_b, new_k, burn);
    MatrixXd res(num_sim, A.cols());
    for(int i = 0; i < num_sim; i++){
        res.row(i) = fr_result.saved_V * steps.row(i).transpose();
    }
    return res; 
}