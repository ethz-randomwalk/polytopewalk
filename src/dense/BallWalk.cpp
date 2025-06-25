#include "BallWalk.hpp"


MatrixXd BallWalk::generateCompleteWalk(const int num_steps, VectorXd& init, const MatrixXd& A, const VectorXd& b, int burn = 0, int seed = -1){
    if (init.rows() != A.cols() || A.rows() != b.rows() ) {
        throw std::invalid_argument("A, b, and init do not match in dimension.");
    }
    VectorXd x = init; 
    int n = x.rows(); 
    int d = A.cols();
    std::mt19937 gen = initializeRNG(seed);
    MatrixXd results = MatrixXd::Zero(num_steps, n);
    int total = (burn + num_steps) * THIN;
    for (int i = 1; i <= total; i++){
        // proposal x_new = x + R /sqrt(d) * Gaussian 
        VectorXd new_x = generateGaussianRVNorm(n, gen) * R/sqrt(d) + x;
        // accept if the proposal is in the polytope
        if (inPolytope(new_x, A, b)){
            x = new_x;
        }
        // if THIN != 1, then record one for every THIN samples 
        if (i % THIN == 0 && i/THIN > burn){
            results.row((int)i/THIN - burn - 1) = x; 
        }
    }
    return results;
}

void BallWalk::printType(){
    cout << "Ball Walk" << endl;
}