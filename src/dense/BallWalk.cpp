#include "BallWalk.hpp"


MatrixXd BallWalk::generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b, int burn = 0){
    int n = x.rows(); 
    int d = A.cols();
    MatrixXd results = MatrixXd::Zero(num_steps, n);
    int total = (burn + num_steps) * THIN;
    for (int i = 1; i <= total; i++){
        VectorXd new_x = generateGaussianRVNorm(n) * R/sqrt(d) + x;
        if (inPolytope(new_x, A, b)){
            x = new_x;
        }
        if (i % THIN == 0 && i/THIN > burn){
            results.row((int)i/THIN - burn - 1) = x; 
        }
    }
    return results;
}

void BallWalk::printType(){
    cout << "Ball Walk" << endl;
}