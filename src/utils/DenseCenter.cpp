#include "DenseCenter.hpp"
#include "DenseLP.hpp"

VectorXd DenseCenter::getInitialPoint(MatrixXd A, VectorXd b){
    int n = A.rows();
    int d = A.cols();
    VectorXd init = VectorXd::Zero(d + 1);
    init(d) = START;

    MatrixXd A2 (n, d + 1);
    for(int i = 0; i < n; i++){
        double norm = A.row(i).norm();
        A2.row(i) << A.row(i)/norm, -1;
        b(i) = b(i)/norm;
    }

    string name = "var_set1";
    Problem lp;
    lp.AddVariableSet(make_shared<ExVariables>(d + 1, name, init));
    lp.AddConstraintSet(make_shared<ExConstraint1>(n, name, A2, b));
    lp.AddCostSet(make_shared<ExCost>(name));
    IpoptSolver ipopt;
    ipopt.SetOption("print_level", 0);
    ipopt.SetOption("sb", "yes");
    ipopt.SetOption("max_iter", MAX_ITER);
    ipopt.SetOption("tol", TOL);
    ipopt.SetOption("s_max", S_MAX);

    ipopt.Solve(lp);
    
    VectorXd sol = lp.GetOptVariables()->GetValues();
    VectorXd temp = sol.head(sol.rows() - 1);
    sol = temp; 
    return sol; 

}