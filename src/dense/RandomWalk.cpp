#include "RandomWalk.hpp"

bool RandomWalk::inPolytope(const VectorXd& vec, const MatrixXd& A, const VectorXd&b){
    return ((A * vec) - b).maxCoeff() <= 0;
}

VectorXd RandomWalk::generateGaussianRVNorm(const int d, std::mt19937& gen){
    VectorXd v(d);
    normal_distribution<double> dis(0.0, 1.0);
    for(int i = 0; i < d; i++){
        v(i) = dis(gen);
    }
    return v/v.norm();
}

MatrixXd RandomWalk::generateCompleteWalk(const int niter, VectorXd& init, 
    const MatrixXd& A, const VectorXd& b, int burnin = 0, int thin = 1, int seed = -1){
    cout << "oops" << endl;
    return MatrixXd::Zero(1,1);
}

void RandomWalk::printType(){
    cout << "oops" << endl;
}

std::mt19937 RandomWalk::initializeRNG(int seed) {
    if (seed != -1) {
        return std::mt19937(seed);
    } else {
        std::random_device rd;
        return std::mt19937(rd());
    }
}