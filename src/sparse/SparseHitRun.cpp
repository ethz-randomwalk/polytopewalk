#include "SparseHitRun.hpp"

double SparseHitAndRun::binarySearch(
    VectorXd direction, 
    VectorXd& x, 
    int k
){
    VectorXd farth = x + R * direction;
    double dist = 0; 

    const int MAXITER = 10000; 
    int iter = 0;

    while(iter < MAXITER){
        dist = (x - farth).norm();
        farth = x + 2 * dist * direction; 
        if (!inPolytope(farth, k)){
            break; 
        }
        iter++;
    }

    if (iter == MAXITER){
        return 0.0;
    }
    VectorXd left = x;
    VectorXd right = farth;
    VectorXd mid = (x + farth)/2;
    while ((left - right).norm() > ERR || !inPolytope(mid, k)){
        mid = (left + right)/2; 
        if (inPolytope(mid, k)){
            left = mid; 
        } else {
            right = mid; 
        }

    }
    return (mid - x).norm();
}

MatrixXd SparseHitAndRun::generateCompleteWalk(
    const int num_steps, 
    const VectorXd& init, 
    const SparseMatrixXd& A, 
    const VectorXd& b, 
    int k,
    int burn = 0
){

    MatrixXd results = MatrixXd::Zero(num_steps, A.cols());
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    SparseLU <SparseMatrixXd> A_solver (A * A.transpose());
    VectorXd x = init; 
    int total = (burn + num_steps) * THIN; 
    for (int i = 1; i <= total; i++){
        VectorXd rand = generateGaussianRV(A.cols());
        VectorXd z = A * rand; 
        z = rand - A.transpose() * A_solver.solve(z);
        z /= z.norm(); 
        double pos_side = binarySearch(z, x, k);
        double neg_side = -binarySearch(-z, x, k);
        double val = dis(gen);
        double random_point = val * (pos_side - neg_side) + neg_side; 
        x = random_point * z + x; 

        if (i % THIN == 0 && i/THIN > burn){
            results.row((int)i/THIN - burn - 1) = x; 
        }
    }
    return results; 

}