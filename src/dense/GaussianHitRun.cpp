#include "GaussianHitRun.hpp"

double GaussianHitAndRun::distance(VectorXd& x, VectorXd&y){
    return (x - y).norm();
}

double GaussianHitAndRun::binarySearch(VectorXd direction, VectorXd& x, const MatrixXd& A, const VectorXd& b){

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


double GaussianHitAndRun::gaussianLogPDF(VectorXd x){

    return -0.5 * (x - MU).transpose() * cov_inv * (x - MU);
}


double GaussianHitAndRun::minF(const VectorXd& v, const VectorXd& x, double l, double u){

    double t_star = - v.transpose() * cov_inv * (x - MU);
    t_star /= (v.transpose() * cov_inv * v);

    if (t_star < l) t_star = l;
    if (t_star > u) t_star = u;

    return gaussianLogPDF(x + v * t_star);

}

MatrixXd GaussianHitAndRun::generateCompleteWalk(const int num_steps, VectorXd& x, const MatrixXd& A, const VectorXd& b, int burn = 0){
    int n = x.rows(); 

    MatrixXd results = MatrixXd::Zero(num_steps, n);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    int total = (burn + num_steps) * THIN; 
    for (int i = 1; i <= total; i++){
        VectorXd new_direct = generateGaussianRVNorm(n);
        double pos_side = binarySearch(new_direct, x, A, b);
        double neg_side = binarySearch(new_direct * -1, x, A, b) * -1;
        double val = dis(gen);
        double random_point = val * (pos_side - neg_side) + neg_side; 
        // the next iterate is uniform on the segment passing x

        VectorXd z = random_point * new_direct + x; 
        double rand = dis(gen); 

        double density1 = gaussianLogPDF(z);
        double density2 = minF(new_direct, x, neg_side, pos_side);

        if (exp(density1) > rand * exp(density2)){
            x = z;
        }
        
        if (i % THIN == 0 && i/THIN > burn){
            results.row((int)i/THIN - burn - 1) = x; 
        }
    }
    return results;
}

void GaussianHitAndRun::printType(){
    cout << "GaussianHitAndRun" << endl;
}