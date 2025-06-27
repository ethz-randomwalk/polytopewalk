#include "SparseBarrierWalk.hpp"

VectorXd SparseBarrierWalk::generateWeight(
    const VectorXd& x, 
    const SparseMatrixXd& A,
    int k
){
    return VectorXd::Ones(A.cols());
}

void SparseBarrierWalk::setDistTerm(int d, int n){
    DIST_TERM = 0;
}

SparseMatrixXd SparseBarrierWalk::generateSlackInverse(const VectorXd& x, int k){
    VectorXd slack_inv (x.rows());
    for(int i = x.rows() - k; i < x.rows(); i++) slack_inv(i) = 1/x(i);
    
    return SparseMatrixXd(slack_inv.asDiagonal());
}

VectorXd SparseBarrierWalk::generateSample(
    const VectorXd& x, 
    const SparseMatrixXd& A, 
    int k,
    std::mt19937& gen
){
    if (k <= 0 || k > A.cols()) {
        throw std::invalid_argument("Parameter k must be between 1 and the number of columns in A.");
    }
    SparseMatrixXd slack_inv = generateSlackInverse(x, k);
    SparseMatrixXd W = SparseMatrixXd(generateWeight(x, A, k).asDiagonal());
    SparseMatrixXd G = slack_inv * W * slack_inv;
    for(int i = 0; i < x.rows() - k; i++) G.coeffRef(i, i) = ERR; 

    SparseMatrixXd G_inv_sqrt = SparseMatrixXd(VectorXd(G.diagonal()).cwiseInverse().cwiseSqrt().asDiagonal());
    
    SparseMatrixXd AG_inv_sqrt = A * G_inv_sqrt;

    VectorXd rand = generateGaussianRV(A.cols(), gen);
    SparseMatrixXd res = AG_inv_sqrt * AG_inv_sqrt.transpose();
    SimplicialLLT<SparseMatrixXd> chol;
    chol.analyzePattern(res);
    chol.factorize(res);
    
    VectorXd z = AG_inv_sqrt * rand;

    z = AG_inv_sqrt.transpose() * chol.solve(z);
    z = G_inv_sqrt * (rand - z);
    z = x + sqrt(DIST_TERM) * z;

    return z; 
}

double SparseBarrierWalk::generateProposalDensity(
    const VectorXd& x, 
    const VectorXd& z, 
    const SparseMatrixXd& A, 
    int k
){
    SparseMatrixXd slack_inv = generateSlackInverse(x, k);
    SparseMatrixXd W = SparseMatrixXd(generateWeight(x, A, k).asDiagonal());
    SparseMatrixXd G = slack_inv * W * slack_inv;
    for(int i = 0; i < x.rows() - k; i++) G.coeffRef(i, i) = ERR; 

    SparseMatrixXd G_inv_sqrt = SparseMatrixXd(VectorXd(G.diagonal()).cwiseInverse().cwiseSqrt().asDiagonal());
    SparseMatrixXd AG_inv_sqrt = A * G_inv_sqrt;
    
    // determinant of S^{-1} W S^{-1}
    double det1 = G.diagonal().array().log().sum();

    SimplicialLLT<SparseMatrixXd> d2; 
    SparseMatrixXd mat = AG_inv_sqrt * AG_inv_sqrt.transpose();
    d2.analyzePattern(mat);
    d2.factorize(mat);

    double det2 = 2 * SparseMatrixXd(d2.matrixL()).diagonal().array().log().sum();
    // -logdet of the matrix g^{-1/2} A^T (A g A^T )^{-1} A g^{-1/2}
    // equals to logdet(g) + logdet(A g A^T) - \logdet(AA^T)
    // but  - \logdet(AA^T) is shared at x or z, so ignored
    double det = det1 + det2; 

    VectorXd diff = z - x;
    VectorXd Qx = A * diff; 
    Qx = A_solver.solve(Qx);
    Qx = A.transpose() * Qx;
    Qx = diff - Qx; 

    double dist = Qx.transpose() * (G * Qx);
    // return the log proposal density
    return 0.5 * det - 0.5/DIST_TERM * dist;
}

MatrixXd SparseBarrierWalk::generateCompleteWalk(
    const int niter, 
    const VectorXd& init, 
    const SparseMatrixXd& A, 
    const VectorXd& b, 
    int k,
    int burnin = 0,
    int thin = 1,
    int seed = -1
){
    if (k < 0 || k > A.cols()) {
        throw std::invalid_argument("Parameter k must be between 0 and the number of columns in A.");
    }

    if (init.rows() != A.cols() || A.rows() != b.rows() ) {
        throw std::invalid_argument("A, b, and init do not match in dimension.");
    }
    int total_samples = (niter - burnin)/thin;
    MatrixXd results = MatrixXd::Zero(total_samples, A.cols());
    std::mt19937 gen = initializeRNG(seed);
    uniform_real_distribution<> dis(0.0, 1.0);

    setDistTerm(A.cols() - A.rows(), k);
    VectorXd x = init;
    A_solver.compute(A * A.transpose());
    for(int i = 1; i <= niter; i++){
        VectorXd z = generateSample(x, A, k, gen);
        if (inPolytope(z, k)){
            double g_x_z = generateProposalDensity(x, z, A, k);
            double g_z_x = generateProposalDensity(z, x, A, k);
            double alpha = min(1.0, exp(g_z_x - g_x_z));
            double val = dis(gen);
            x = val < alpha ? z : x; 
        }
        if (i > burnin && (i - burnin) % thin == 0){
            results.row((int)((i - burnin)/thin - 1)) = x; 
        }
    }

    return results; 

}