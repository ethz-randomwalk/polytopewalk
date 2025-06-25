#include "DikinWalk.hpp"

void DikinWalk::setDistTerm(int d, int n){
    DIST_TERM = R*R/d;
}

VectorXd DikinWalk::generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    int d = b.rows();
    return VectorXd::Ones(d);
}

void DikinWalk::printType(){
    cout << "Dikin Walk" << endl;
}