#include "BallWalk.hpp"


MatrixXd BallWalk::generateCompleteWalk(const int niter, VectorXd& init, const MatrixXd& A, const VectorXd& b, 
    int burnin = 0, int thin = 1, int seed = -1){

    if (init.rows() != A.cols() || A.rows() != b.rows() ) {
        throw std::invalid_argument("A, b, and init do not match in dimension.");
    }
    VectorXd x = init; 
    int n = x.rows(); 
    int d = A.cols();
    std::mt19937 gen = initializeRNG(seed);
    int total_samples = (niter - burnin)/thin;
    MatrixXd results = MatrixXd::Zero(total_samples, n);
    for (int i = 1; i <= niter; i++){
        // proposal x_new = x + R /sqrt(d) * Gaussian 
        VectorXd new_x = generateGaussianRVNorm(n, gen) * R/sqrt(d) + x;
        // accept if the proposal is in the polytope
        if (inPolytope(new_x, A, b)){
            x = new_x;
        }
        // if thin != 1, then record one for every thin samples 
        if (i > burnin && (i - burnin) % thin == 0){
            results.row((int)((i - burnin)/thin - 1)) = x; 
        }
    }
    return results;
}

void BallWalk::printType(){
    cout << "Ball Walk" << endl;
}