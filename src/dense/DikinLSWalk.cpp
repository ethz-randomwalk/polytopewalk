#include "DikinLSWalk.hpp"

void DikinLSWalk::setDistTerm(int d, int n){
    w_i = VectorXd::Ones(n); 
    double q = 2.0 * (1.0 + log(n));
    double term = (1.0 + q) * (1.0 + q * q);
    DIST_TERM = R*R/term; 
}

void DikinLSWalk::generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b){

    double q = 2.0 * (1.0 + log(A.rows()));
    double alpha = 1.0 - (2.0/q);

    if (w_i.coeffRef(0) == -1 || w_i.rows() != A.rows()){
        w_i = VectorXd::Ones(A.rows());
    }

    generateSlack(x, A, b);
    DiagonalMatrix<double, Dynamic> slack_inv = slack.cwiseInverse().asDiagonal();
    MatrixXd A_x = slack_inv * A; 

    DiagonalMatrix<double, Dynamic> W;
    MatrixXd WAX (A.rows(), A.cols());
    VectorXd term2a (A.rows());
    VectorXd term2b (A.rows());
    VectorXd term2(A.rows());
    VectorXd gradient (A.rows()); 
    VectorXd proposal (A.rows()); 
    VectorXd term3 (A.rows());
    VectorXd error = 0.00001 * VectorXd::Ones(A.rows());

    // gradient descent to compute LS weights
    for(int i = 0; i < MAXITER; i++){
        W = VectorXd(w_i.array().pow(alpha * 0.5)).asDiagonal();
        term2a = alpha * w_i.cwiseInverse();

        WAX = W * A_x;
        // leverage score based on previous W
        term2b = (WAX * (WAX.transpose() * WAX).inverse()).cwiseProduct(WAX).rowwise().sum();

        term2 = term2a.cwiseProduct(term2b);
        
        gradient =  (alpha) * VectorXd::Ones(A.rows()) - term2;
        if(gradient.norm() < GRADLIM){
            break;
        }
        w_i = (w_i - STEPSIZE * gradient).cwiseMax(error);
    }
    weights = w_i.asDiagonal();
    
}

void DikinLSWalk::printType(){
    cout << "DikinLSWalk" << endl;
}