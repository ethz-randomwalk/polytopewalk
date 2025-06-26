#include "HitRun.hpp"

double HitAndRun::distance(VectorXd& x, VectorXd&y){
    return (x - y).norm();
}

double HitAndRun::binarySearch(VectorXd direction, VectorXd& x, const MatrixXd& A, const VectorXd& b){

    VectorXd farth = x + R * direction;
    double dist = 0; 

    while(true){
        dist = distance(x, farth);
        farth = x + 2 * dist * direction; 
        if (!inPolytope(farth, A, b)){
            break; 
        }
    }
    VectorXd left = x;
    VectorXd right = farth;
    VectorXd mid = (x + farth)/2;

    while (distance(left, right) > ERR || ! inPolytope(mid, A, b)){
        mid = (left + right)/2; 
        if (inPolytope(mid, A, b)){
            left = mid; 
        } else {
            right = mid; 
        }
    }
    // return the distance bewteen the intersection of direction and polytope
    // and x
    return distance(mid, x);
}

MatrixXd HitAndRun::generateCompleteWalk(const int niter, VectorXd& init, const MatrixXd& A, const VectorXd& b, int burnin = 0, int thin = 1, int seed = -1){
    if (init.rows() != A.cols() || A.rows() != b.rows() ) {
        throw std::invalid_argument("A, b, and init do not match in dimension.");
    }

    VectorXd x = init; 
    
    int n = x.rows(); 
    int total_samples = (niter - burnin)/thin;
    MatrixXd results = MatrixXd::Zero(total_samples, n);
    std::mt19937 gen = initializeRNG(seed);
    uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 1; i <= niter; i++){
        VectorXd new_direct = generateGaussianRVNorm(n, gen);
        double pos_side = binarySearch(new_direct, x, A, b);
        double neg_side = binarySearch(new_direct * -1, x, A, b) * -1;
        double val = dis(gen);
        double random_point = val * (pos_side - neg_side) + neg_side; 
        // the next iterate is uniform on the segment passing x
        x = random_point * new_direct + x; 
        
        if (i > burnin && (i - burnin) % thin == 0){
            results.row((int)((i - burnin)/thin - 1)) = x; 
        }
    }
    return results;
}

void HitAndRun::printType(){
    cout << "HitAndRunWalk" << endl;
}