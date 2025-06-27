#include "SparseBallWalk.hpp"

MatrixXd SparseBallWalk::generateCompleteWalk(
    const int niter, 
    const VectorXd& init, 
    const SparseMatrixXd& A, 
    const VectorXd& b, 
    int k, 
    int burnin = 0,
    int thin = 1,
    int seed = -1
){
    if (k <= 0 || k > A.cols()) {
        throw std::invalid_argument("Parameter k must be between 1 and the number of columns in A.");
    }
    if (init.rows() != A.cols() || A.rows() != b.rows() ) {
        throw std::invalid_argument("A, b, and init do not match in dimension.");
    }
    int total_samples = (niter - burnin)/thin;
    MatrixXd results = MatrixXd::Zero(total_samples, A.cols());

    SparseLU<SparseMatrixXd> A_solver (A * A.transpose());
    SparseMatrixXd I = SparseMatrixXd(VectorXd::Ones(A.cols()).asDiagonal());

    std::mt19937 gen = initializeRNG(seed);

    VectorXd x = init;
    int d = A.cols() - A.rows();
    for (int i = 1; i <= niter; i++){
        VectorXd rand = generateGaussianRV(A.cols(), gen); 
        VectorXd z;
        z = A * rand; 
        z = rand - A.transpose() * A_solver.solve(z);
        z /= z.norm(); 
        z = R/sqrt(d) * z + x; 

        if (inPolytope(z, k)){
            x = z;
        } 
        if (i > burnin && (i - burnin) % thin == 0){
            results.row((int)((i - burnin)/thin - 1)) = x; 
        }
    }
    return results; 
}