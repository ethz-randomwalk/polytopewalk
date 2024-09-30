#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/FullWalkRun.hpp"
#include <cstring>

struct sparse_polytope{
    SparseMatrixXd A;
    VectorXd b; 
    int k;
};

sparse_polytope generate_simplex(){
    SparseMatrixXd simplex_A (1, 3);
    simplex_A.coeffRef(0, 0) = 1;
    simplex_A.coeffRef(0, 1) = 1;
    simplex_A.coeffRef(0, 2) = 1;
    VectorXd simplex_b (1);
    simplex_b << 1;
    sparse_polytope result; 
    result.A = simplex_A;
    result.b = simplex_b;
    result.k = 3; 
    return result;
}

sparse_polytope generate_hc(){
    SparseMatrixXd hc_A (4, 6);
    hc_A.coeffRef(0, 0) = 1;
    hc_A.coeffRef(0, 2) = 1;
    hc_A.coeffRef(1, 1) = 1;
    hc_A.coeffRef(1, 3) = 1;
    hc_A.coeffRef(2, 0) = -1;
    hc_A.coeffRef(2, 4) = 1;
    hc_A.coeffRef(3, 1) = -1;
    hc_A.coeffRef(3, 5) = 1; 

    VectorXd hc_b (4);
    hc_b << 1, 1, 1, 1;
    sparse_polytope result; 
    result.A = hc_A;
    result.b = hc_b;
    result.k = 4; 
    return result;
}

sparse_polytope generate_birkhoff(){
    SparseMatrixXd birk_A (3, 4);
    birk_A.coeffRef(0, 0) = 1;
    birk_A.coeffRef(0, 1) = 1;
    birk_A.coeffRef(1, 2) = 1;
    birk_A.coeffRef(1, 3) = 1;
    birk_A.coeffRef(2, 0) = 1;
    birk_A.coeffRef(2, 2) = 1; 

    VectorXd birk_b (3);
    birk_b << 1, 1, 1;
    sparse_polytope result; 
    result.A = birk_A;
    result.b = birk_b;
    result.k = 4; 
    return result;
}

sparse_polytope simplex = generate_simplex();
sparse_polytope hc = generate_hc();
sparse_polytope birk = generate_birkhoff();

TEST_CASE( "Check Facial Reduction Algorithm", "[require]" ) {
    
    FacialReduction fr;

    res simplex_dense = fr.reduce(simplex.A, simplex.b, simplex.k, false);
    res hc_dense = fr.reduce(hc.A, hc.b, hc.k, false);
    res birk_dense = fr.reduce(birk.A, birk.b, birk.k, false);

    int dense_A_row, dense_A_col, dense_b_row;
    dense_A_row = simplex_dense.dense_A.rows();
    dense_A_col = simplex_dense.dense_A.cols();
    dense_b_row = simplex_dense.dense_b.rows();

    REQUIRE(((dense_A_row == 3) && (dense_A_col == 2)));
    REQUIRE(dense_b_row == 3);

    dense_A_row = hc_dense.dense_A.rows();
    dense_A_col = hc_dense.dense_A.cols();
    dense_b_row = hc_dense.dense_b.rows();

    REQUIRE(((dense_A_row == 4) && (dense_A_col == 2)));
    REQUIRE(dense_b_row == 4);

    dense_A_row = birk_dense.dense_A.rows();
    dense_A_col = birk_dense.dense_A.cols();
    dense_b_row = birk_dense.dense_b.rows();

    REQUIRE(((dense_A_row == 4) && (dense_A_col == 1)));
    REQUIRE(dense_b_row == 4);

    MatrixXd A1 (6, 3);
    A1 << 1, 1, 0, -1, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
    MatrixXd temp (6, 9);
    temp << A1, VectorXd::Ones(6).asDiagonal().toDenseMatrix();
    SparseMatrixXd SA1 = temp.sparseView();
    VectorXd b1(6);
    b1 << 1, -1, 1, 1, 1, 1;
    
    
    res test1a = fr.reduce(SA1, b1, 6, true);
    REQUIRE((test1a.sparse_A.rows() == 5 && test1a.sparse_A.cols() == 7));
    res test1b = fr.reduce(SA1, b1, 6, false);
    dense_A_row = test1b.dense_A.rows();
    dense_A_col = test1b.dense_A.cols();
    REQUIRE((dense_A_row == 4 && dense_A_col == 2));
    MatrixXd A2(6,3);
    A2 << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
    MatrixXd temp2 (6, 9);
    temp2 << A2, VectorXd::Ones(6).asDiagonal().toDenseMatrix();
    SparseMatrixXd SA2 = temp2.sparseView();
    VectorXd b2(6);
    b2 << 1, 1, 0, 0, 0, 0;

    res test2a = fr.reduce(SA2, b2, 6, true);
    int sparse_rows = test2a.sparse_A.rows();
    int sparse_cols = test2a.sparse_A.cols();

    REQUIRE_THAT(sparse_rows, Catch::Matchers::WithinAbs(4, 0.01));
    REQUIRE_THAT(sparse_cols, Catch::Matchers::WithinAbs(5, 0.01));
    res test2b = fr.reduce(SA2, b2, 6, false);
    dense_A_row = test2b.dense_A.rows();
    dense_A_col = test2b.dense_A.cols();
    REQUIRE((dense_A_row == 2 && dense_A_col == 1));


    MatrixXd A3(4,2);
    A3 << 1, 0, -1, 0, 0, 1, 0, -1;
    MatrixXd temp3 (4, 6);
    temp3 << A3, VectorXd::Ones(4).asDiagonal().toDenseMatrix();
    SparseMatrixXd SA3 = temp3.sparseView();
    VectorXd b3(4);
    b3 << 1, 0, 1, 0;

    res test3a = fr.reduce(SA3, b3, 4, true);
    int sparse_A_rows = test3a.sparse_A.rows(); 
    int sparse_A_cols = test3a.sparse_A.cols(); 

    REQUIRE_THAT(sparse_A_rows, Catch::Matchers::WithinAbs(4, 0.0001));
    REQUIRE_THAT(sparse_A_cols, Catch::Matchers::WithinAbs(6, 0.0001));
    res test3b = fr.reduce(SA3, b3, 4, false);
    dense_A_row = test3b.dense_A.rows();
    dense_A_col = test3b.dense_A.cols();
    REQUIRE((dense_A_row == 4 && dense_A_col == 2));


}
