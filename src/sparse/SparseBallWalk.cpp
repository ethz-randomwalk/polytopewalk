#include "SparseBallWalk.hpp"

MatrixXd SparseBallWalk::generateCompleteWalk(
    const int num_steps, 
    const VectorXd& init, 
    const SparseMatrixXd& A, 
    const VectorXd& b, 
    int k, 
    int burn = 0,
    int seed = -1
){
    if (k < 0 || k > A.cols()) {
        throw std::invalid_argument("Parameter k must be between 0 and the number of columns in A.");
    }
    if (init.rows() != A.cols() || A.rows() != b.rows() ) {
        throw std::invalid_argument("A, b, and init do not match in dimension.");
    }
    
    MatrixXd results = MatrixXd::Zero(num_steps, A.cols());

    SparseLU<SparseMatrixXd> A_solver (A * A.transpose());
    SparseMatrixXd I = SparseMatrixXd(VectorXd::Ones(A.cols()).asDiagonal());

    std::mt19937 gen = initializeRNG(seed);

    VectorXd x = init;
    int d = A.cols() - A.rows();
    int total = (burn + num_steps) * THIN;
    for (int i = 1; i <= total; i++){
        VectorXd rand = generateGaussianRV(A.cols(), gen); 
        VectorXd z;
        z = A * rand; 
        z = rand - A.transpose() * A_solver.solve(z);
        z /= z.norm(); 
        z = R/sqrt(d) * z + x; 

        if (inPolytope(z, k)){
            x = z;
        } 
        if (i % THIN == 0 && i/THIN > burn){
            results.row((int)i/THIN - burn - 1) = x; 
        }
    }
    return results; 
}