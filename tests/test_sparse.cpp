
#include "utils/FullWalkRun.hpp"

int main(){
    SparseMatrixXd A (3, 4);
    A.coeffRef(0, 0) = 1;
    A.coeffRef(0, 1) = 1;
    A.coeffRef(1, 2) = 1;
    A.coeffRef(1, 3) = 1;
    A.coeffRef(2, 0) = 1;
    A.coeffRef(2, 2) = 1; 

    VectorXd b (3);
    b << 1, 1, 1;
    
    VectorXd x(4);
    x << 0.5, 0.5, 0.5, 0.5;

    SparseMatrixXd C (3, 3);
    C.coeffRef(0, 0) = 2;
    C.coeffRef(0, 1) = 1;
    C.coeffRef(0, 2) = 0.5;
    C.coeffRef(1, 0) = 3;
    C.coeffRef(1, 1) = 2.5;
    C.coeffRef(1, 2) = 4;
    C.coeffRef(2, 0) = 1.5;
    C.coeffRef(2, 1) = 1;
    C.coeffRef(2, 2) = 1;


    SparseVaidyaWalk vaidya (0.6);
    MatrixXd res = vaidya.generateCompleteWalk(10, x, A, b, 3, 0);
}