#include "DikinWalk.hpp"

void DikinWalk::setDistTerm(int d, int n){
    DIST_TERM = R*R/d;
}

void DikinWalk::generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    int d = b.rows();
    weights = VectorXd::Ones(d).asDiagonal();
}

void DikinWalk::printType(){
    cout << "Dikin Walk" << endl;
}