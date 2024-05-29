#include "FacialReduction.hpp"


z_res FacialReduction::findZ(const SparseMatrixXd& A, const VectorXd& b, int x_dim){
    int n = A.rows();
    int d = A.cols();

    z_res ans;
    ans.found_sol = false;

    SparseMatrixXd ineqA(n + 1, d-x_dim);
    for(int i = x_dim; i < d; i++){
        ineqA.col(i - x_dim) = -1 * A.col(i);
    }
    ineqA = ineqA.transpose();
    ineqA.col(n) = (-1 * VectorXd::Ones(d-x_dim)).sparseView();

    SparseMatrixXd eqA(n + 1, x_dim + 2);
    for(int j = 0; j < x_dim; j++){
        eqA.col(j) = A.col(j);
    }
    eqA.col(x_dim + 1) = b.sparseView();
    VectorXd eqb = VectorXd::Zero(eqA.cols());

    for(int i = global_index; i < d; i++){
        eqA.col(x_dim) = A.col(i);
        eqb(x_dim) = 1.0;
        eqA = eqA.transpose();

        SparseQR<SparseMatrixXd, COLAMDOrdering<SparseMatrix<double>::StorageIndex>> solver (eqA.block(0, 0, x_dim + 2, n));
        VectorXd init = solver.solve(eqb);
        VectorXd temp_solve (init.rows() + 1);
        temp_solve << init, 0; 
        double delta = (ineqA * temp_solve).maxCoeff();
        VectorXd temp (init.rows() + 1);
        temp << init, delta;
        init = temp; 

        if((eqA * init - eqb).cwiseAbs().maxCoeff() > ERR_LP){
            eqA = eqA.transpose();
            global_index ++; 
            continue;
        }

        if(init(init.rows() - 1) <= ERR_LP){
            VectorXd temp = init.head(init.rows() - 1);
            init = temp;

            ans.found_sol = true; 
            ans.z = (A.transpose() * init);
            return ans;
        }

        string name = "var_set1";
        Problem lp;
        lp.AddVariableSet(make_shared<SparseExVariables>(n + 1, name, init));
        lp.AddConstraintSet(make_shared<SparseExConstraint1>(ineqA.rows(),name,ineqA, ERR_LP));
        lp.AddConstraintSet(make_shared<SparseExConstraint2>(eqA.rows(),name,eqA,eqb, ERR_LP));
        lp.AddCostSet(make_shared<SparseExCost>(name));
        IpoptSolver ipopt;
        ipopt.SetOption("print_level", 0);
        ipopt.SetOption("sb", "yes");
        ipopt.SetOption("max_iter", MAX_ITER);
        ipopt.SetOption("tol", TOL);
        ipopt.SetOption("s_max", S_MAX);
        ipopt.Solve(lp);

        VectorXd sol = lp.GetOptVariables()->GetValues();

        if (ipopt.GetReturnStatus() != 0){
            global_index ++;
            eqA = eqA.transpose();
            continue; 
        }

        if (sol(sol.rows() - 1) <= ERR_LP){
            VectorXd temp = sol.head(sol.rows() - 1);
            sol = temp;
            ans.found_sol = true; 
            ans.z = (A.transpose() * sol);
            return ans;
        }
        global_index ++;
        eqA = eqA.transpose();

    }
    return ans;
}

SparseMatrixXd FacialReduction::pickV(const VectorXd& z, int x_dim){
    int d = z.rows();
    vector<T> indices;
    for(int i = 0; i < x_dim; i++){
        indices.push_back(T(indices.size(), i, 1)); 
    }
    for(int i = x_dim; i < d; i++){
         if(z(i) < ERR_DC) indices.push_back(T(indices.size(), i, 1)); 
    }


    SparseMatrixXd mat(indices.size(), d);
    mat.setFromTriplets(indices.begin(), indices.end());
    return mat.transpose();
}

SparseMatrixXd FacialReduction::pickP(const SparseMatrixXd& AV){
    SparseQR<SparseMatrixXd, NaturalOrdering<SparseMatrix<double>::StorageIndex>> solver;
    solver.compute(AV.transpose());
    SparseMatrixXd R = solver.matrixR();

    vector<T> indices;
    for (int i = 0; i < min(R.cols(), R.rows()); i++){
        if (abs(R.coeffRef(i, i)) > ERR_DC){
            indices.push_back(T(indices.size(), solver.colsPermutation().indices()(i), 1));
        }
    }
    SparseMatrixXd proj (indices.size(), AV.rows());
    proj.setFromTriplets(indices.begin(), indices.end());
    return proj; 
}

fr_res FacialReduction::entireFacialReductionStep(SparseMatrixXd A, VectorXd b, int x_dim){
    z_res z_ans = findZ(A, b, x_dim);

    if(!z_ans.found_sol){
        fr_res ans;
        ans.A = A;
        ans.b = b; 
        return ans; 
    }
    SparseMatrixXd V = pickV(z_ans.z, x_dim);
    savedV = savedV * V; 
    SparseMatrixXd AV = A * V;
    SparseMatrixXd P = pickP(AV);
    A = P * AV;
    b = P * b; 
    return entireFacialReductionStep(A, b, x_dim);
}

res FacialReduction::reduce(SparseMatrixXd A, VectorXd b, int k, bool sparse){
    int x_dim = A.cols() - k; 
    savedV = SparseMatrixXd(VectorXd::Ones(A.cols()).asDiagonal());
    global_index = x_dim; 
    //remove dependent rows
    SparseMatrixXd P = pickP(A);
    A = P * A; 
    b = P * b; 
    fr_res result = entireFacialReductionStep(A, b, x_dim);
    res final_res; 

    final_res.sparse_A = result.A;
    final_res.sparse_b = result.b;
    final_res.saved_V = savedV;

    if(!sparse){
        HouseholderQR <MatrixXd> qr(result.A.cols(), result.A.rows());
        qr.compute(MatrixXd(result.A.transpose()));
        MatrixXd Q = qr.householderQ();
        MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
        int d = R.rows();
        int n = R.cols();

        MatrixXd newR = R.block(0, 0, R.cols(), R.cols());
        VectorXd z1 = newR.transpose().inverse() * result.b;

        MatrixXd Q1 = Q.block(0, 0, Q.rows(), n);
        MatrixXd Q2 = Q.block(0, n, Q.rows(), d - n);
        MatrixXd reduced_A = -1 * Q2.block(x_dim, 0, d - x_dim, d - n);
        VectorXd reduced_b = (Q1 * z1).tail(d - x_dim);

        final_res.dense_A = reduced_A;
        final_res.dense_b = reduced_b;
        
        final_res.z1 = z1;
        final_res.Q = Q;
    }
    return final_res;
}

