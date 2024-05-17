#include "JohnWalk.hpp"

void JohnWalk::setDistTerm(int d, int n){
    w_i = VectorXd::Ones(n); 
    DIST_TERM = R*R/(pow(d, 1.5));
}

void JohnWalk::generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b ){
    double alpha = 1 - 1/(log2(2 * A.rows() / A.cols()));
    double beta = (double)A.cols() / (2 * A.rows());

    generateSlack(x, A, b);
    DiagonalMatrix<double, Dynamic> slack_inv = slack.cwiseInverse().asDiagonal();

    if (w_i.coeffRef(0) == -1 || w_i.rows() != A.rows()){
        w_i = VectorXd::Ones(A.rows());
    }

    MatrixXd A_x = slack_inv * A; 

    DiagonalMatrix<double, Dynamic> W;
    MatrixXd WAX (A.rows(), A.cols());
    VectorXd gradient (A.rows());
    VectorXd score; 

    VectorXd beta_ones = beta * VectorXd::Ones(A.rows());
    VectorXd next_weight = w_i; 

    for(int i = 0; i < MAXITER; i++){
        w_i = next_weight; 

        W = VectorXd(w_i.array().pow(alpha * 0.5)).asDiagonal();
        WAX = W * A_x;
        score = (WAX * (WAX.transpose() * WAX).inverse()).cwiseProduct(WAX).rowwise().sum();

        next_weight = 0.5 * (w_i + score + beta_ones).cwiseMax(beta_ones);
        if((next_weight - w_i).cwiseAbs().maxCoeff() < LIM){
            break;
        }
    }

    weights = w_i.asDiagonal();

}


void JohnWalk::printType(){
    cout << "John Walk" << endl;
}