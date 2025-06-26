#include "SparseRandomWalk.hpp"

VectorXd SparseRandomWalk::generateGaussianRV(int d, std::mt19937& gen){
    normal_distribution<double> dis(0.0, 1.0);
    VectorXd v(d);
    for(int i = 0; i < d; i++){
        v(i) = dis(gen);
    }
    return v;
}

MatrixXd SparseRandomWalk::generateCompleteWalk(
    const int niter,
    const VectorXd& init, 
    const SparseMatrixXd& A,
    const VectorXd& b, 
    int k,
    int burnin = 0,
    int thin = 1,
    int seed = -1
){
    cout << "Oops" << endl;
    return MatrixXd::Zero(1,1);

}

bool SparseRandomWalk::inPolytope(
    const VectorXd&z, 
    int k
){
    return z.tail(k).minCoeff() >= 0; 
}


std::mt19937 SparseRandomWalk::initializeRNG(int seed) {
    if (seed != -1) {
        return std::mt19937(seed);
    } else {
        std::random_device rd;
        return std::mt19937(rd());
    }
}

