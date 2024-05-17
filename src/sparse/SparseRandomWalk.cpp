#include "SparseRandomWalk.hpp"

VectorXd SparseRandomWalk::generateGaussianRV(int d){
    VectorXd v(d);
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dis(0.0, 1.0);
    for(int i = 0; i < d; i++){
        v(i) = dis(gen);
    }
    return v;
}

MatrixXd SparseRandomWalk::generateCompleteWalk(
    const int num_steps,
    const VectorXd& init, 
    const SparseMatrixXd& A,
    const VectorXd& b, 
    int k,
    int burn = 0
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