#include "SparseGaussianHitRun.hpp"

double SparseGaussianHitAndRun::binarySearch(
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

double SparseGaussianHitAndRun::gaussianLogPDF(VectorXd x){
    return -0.5 * (x-MU).transpose() *  chol.solve(x - MU);
}


double SparseGaussianHitAndRun::minF(VectorXd& v, VectorXd& x, double l, double u){
    double t_star = - v.transpose() * chol.solve(x - MU);
    t_star /= (v.transpose() * chol.solve(v));

    if (t_star < l) t_star = l;
    if (t_star > u) t_star = u;

    return gaussianLogPDF(x + v * t_star);
}

MatrixXd SparseGaussianHitAndRun::generateCompleteWalk(
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
        VectorXd new_direct = generateGaussianRV(A.cols());
        VectorXd new_direct_proj = A*new_direct; 
        new_direct_proj = new_direct - A.transpose() * A_solver.solve(new_direct_proj);
        new_direct_proj /= new_direct_proj.norm(); 
        double pos_side = binarySearch(new_direct_proj, x, k);
        double neg_side = -binarySearch(-new_direct_proj, x, k);
        double val = dis(gen);
        double random_point = val * (pos_side - neg_side) + neg_side; 
        
        VectorXd z = random_point * new_direct_proj + x; 
        double rand = dis(gen); 

        double density1 = gaussianLogPDF(z);
        double density2 = minF(new_direct, x, neg_side, pos_side);

        if (exp(density1) > rand *  exp(density2)){
            x = z;
        }

        if (i % THIN == 0 && i/THIN > burn){
            results.row((int)i/THIN - burn - 1) = x; 
        }
    }
    return results; 

}