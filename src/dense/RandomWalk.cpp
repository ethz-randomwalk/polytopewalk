#include "RandomWalk.hpp"

bool RandomWalk::inPolytope(const VectorXd& vec, const MatrixXd& A, const VectorXd&b){
    return ((A * vec) - b).maxCoeff() <= 0;
}

VectorXd RandomWalk::generateGaussianRVNorm(const int d){
    VectorXd v(d);
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dis(0.0, 1.0);
    for(int i = 0; i < d; i++){
        v(i) = dis(gen);
    }
    return v/v.norm();
}

MatrixXd RandomWalk::generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b, int burn = 0){
    cout << "oops" << endl;
    return MatrixXd::Zero(1,1);
}

void RandomWalk::printType(){
    cout << "oops" << endl;
}