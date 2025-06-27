#include "BarrierWalk.hpp"

void BarrierWalk::setDistTerm(int d, int n){
    DIST_TERM = R*R/n;
}

VectorXd BarrierWalk::generateGaussianRV(int d, std::mt19937& gen){
    VectorXd v(d);
    normal_distribution<double> dis(0.0, 1.0);
    for(int i = 0; i < d; i++){
        v(i) = dis(gen);
    }
    return v;
}

VectorXd BarrierWalk::generateSlack(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    return (b - (A * x));
}

double BarrierWalk::localNorm(VectorXd v, const MatrixXd& m){
    return ((v.transpose() * m) * v)(0);
}

VectorXd BarrierWalk::generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    // always overwrite
    return VectorXd::Ones(1); 
}

MatrixXd BarrierWalk::generateHessian(const VectorXd& x, const MatrixXd& A, const VectorXd& b){
    VectorXd weights = generateWeight(x, A, b);
    VectorXd slack = generateSlack(x, A, b);
    VectorXd slack_inv = slack.cwiseInverse();
    DiagonalMatrix<double, Dynamic> middle = slack_inv.cwiseProduct(weights).cwiseProduct(slack_inv).asDiagonal();
    MatrixXd hess = A.transpose() * middle * A;
    return hess;
}

VectorXd BarrierWalk::generateSample(const VectorXd& x, const MatrixXd& A, const VectorXd& b, std::mt19937& gen){
    uniform_real_distribution<> dis(0.0, 1.0);

    MatrixXd hess = generateHessian(x, A, b);
    // cholesky decomposition to compute inverse of hess
    LLT<MatrixXd> cholesky1(hess);
    MatrixXd L = cholesky1.matrixL();
    FullPivLU<MatrixXd> lu(L);
    VectorXd direction = generateGaussianRV(x.rows(), gen);
    VectorXd prop = x + sqrt(DIST_TERM) * (lu.solve(direction));

    if(!inPolytope(prop, A, b)){
        return x; 
    }
    
    double det = L.diagonal().array().log().sum(); 
    double dist = -(0.5/DIST_TERM) * localNorm(x - prop, hess);
    double g_x_z = det + dist; 

    hess = generateHessian(prop, A, b);
    LLT<MatrixXd> cholesky2(hess);
    L = cholesky2.matrixL();
    det = L.diagonal().array().log().sum(); 
    dist = -(0.5/DIST_TERM) * localNorm(x - prop, hess);
    double g_z_x = det + dist;  

    // accept reject step
    double alpha = min(1.0, exp(g_z_x-g_x_z));
    double val = dis(gen);
    prop = val < alpha ? prop : x;
    return prop; 
}

MatrixXd BarrierWalk::generateCompleteWalk(const int niter, VectorXd& init, const MatrixXd& A, const VectorXd& b, int burnin = 0, int thin = 1, int seed = -1){
    
    if (init.rows() != A.cols() || A.rows() != b.rows() ) {
        throw std::invalid_argument("A, b, and init do not match in dimension.");
    }
    int total_samples = (niter - burnin)/thin;
    MatrixXd results = MatrixXd::Zero(total_samples, A.cols());
    std::mt19937 gen = initializeRNG(seed);

    VectorXd x = init; 

    setDistTerm(A.cols(), A.rows());
    for(int i = 1; i <= niter; i++){
        VectorXd prop = generateSample(x, A, b, gen);
        x = prop; 

       if (i > burnin && (i - burnin) % thin == 0){
            results.row((int)((i - burnin)/thin - 1)) = x; 
        }
    }
    return results;
}