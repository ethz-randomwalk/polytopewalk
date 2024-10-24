#include "BarrierWalk.hpp"

void BarrierWalk::setDistTerm(int d, int n){
    DIST_TERM = R*R/n;
}

VectorXd BarrierWalk::generateGaussianRV(int d){
    VectorXd v(d);
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dis(0.0, 1.0);
    for(int i = 0; i < d; i++){
        v(i) = dis(gen);
    }
    return v;
}

void BarrierWalk::generateSlack(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    slack = (b - (A * x));
}

double BarrierWalk::localNorm(VectorXd v, const MatrixXd& m){
    return ((v.transpose() * m) * v)(0);
}

void BarrierWalk::generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    // always overwrite
}

void BarrierWalk::generateHessian(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    generateWeight(x, A, b);
    generateSlack(x, A, b);
    VectorXd slack_inv = slack.cwiseInverse();
    DiagonalMatrix<double, Dynamic> middle = slack_inv.cwiseProduct(weights.diagonal()).cwiseProduct(slack_inv).asDiagonal();
    hess = A.transpose() * middle * A;
}

void BarrierWalk::generateSample(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    generateHessian(x, A, b); // sets global hess
    // cholesky decomposition to compute inverse of hess
    LLT<MatrixXd> cholesky1(hess);
    MatrixXd L = cholesky1.matrixL();
    FullPivLU<MatrixXd> lu(L);
    VectorXd direction = generateGaussianRV(x.rows());
    prop = x + sqrt(DIST_TERM) * (lu.solve(direction));

    if(!inPolytope(prop, A, b)){
        prop = x;
        return; 
    }
    double det = L.diagonal().array().log().sum(); 
    double dist = -(0.5/DIST_TERM) * localNorm(x - prop, hess);
    double g_x_z = det + dist; 

    generateHessian(prop, A, b);
    LLT<MatrixXd> cholesky2(hess);
    L = cholesky2.matrixL();
    det = L.diagonal().array().log().sum(); 
    dist = -(0.5/DIST_TERM) * localNorm(x - prop, hess);
    double g_z_x = det + dist;  

    // accept reject step
    double alpha = min(1.0, exp(g_z_x-g_x_z));
    double val = dis(gen);
    prop = val < alpha ? prop : x;
}

MatrixXd BarrierWalk::generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b, int burn = 0){
    MatrixXd results = MatrixXd::Zero(num_steps, A.cols());

    setDistTerm(A.cols(), A.rows());
    int total = (burn + num_steps) * THIN; 
    for(int i = 1; i <= total; i++){
        generateSample(x, A, b);
        x = prop; 

        if (i % THIN == 0 && i/THIN > burn){
            results.row((int)i/THIN - burn - 1) = x.transpose(); 
        }
    }
    return results;
}