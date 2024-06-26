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

    REQUIRE(((simplex_dense.dense_A.rows() == 3) && (simplex_dense.dense_A.cols() == 2)));
    REQUIRE(simplex_dense.dense_b.rows() == 3);
    REQUIRE(((hc_dense.dense_A.rows() == 4) && (hc_dense.dense_A.cols() == 2)));
    REQUIRE(hc_dense.dense_b.rows() == 4);
    REQUIRE(((birk_dense.dense_A.rows() == 4) && (birk_dense.dense_A.cols() == 1)));
    REQUIRE(birk_dense.dense_b.rows() == 4);

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
    REQUIRE((test1b.dense_A.rows() == 4 && test1b.dense_A.cols() == 2));

    MatrixXd A2(6,3);
    A2 << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
    MatrixXd temp2 (6, 9);
    temp2 << A2, VectorXd::Ones(6).asDiagonal().toDenseMatrix();
    SparseMatrixXd SA2 = temp2.sparseView();
    VectorXd b2(6);
    b2 << 1, 1, 0, 0, 0, 0;

    res test2a = fr.reduce(SA2, b2, 6, true);
    REQUIRE((test2a.sparse_A.rows() == 4 && test2a.sparse_A.cols() == 5));
    res test2b = fr.reduce(SA2, b2, 6, false);
    REQUIRE((test2b.dense_A.rows() == 2 && test2b.dense_A.cols() == 1));

    MatrixXd A3(4,2);
    A3 << 1, 0, -1, 0, 0, 1, 0, -1;
    MatrixXd temp3 (4, 6);
    temp3 << A3, VectorXd::Ones(4).asDiagonal().toDenseMatrix();
    SparseMatrixXd SA3 = temp3.sparseView();
    VectorXd b3(4);
    b3 << 1, 0, 1, 0;

    res test3a = fr.reduce(SA3, b3, 4, true);
    REQUIRE((test3a.sparse_A.rows() == 4 && test3a.sparse_A.cols() == 6));
    res test3b = fr.reduce(SA3, b3, 4, false);
    REQUIRE((test3b.dense_A.rows() == 4 && test3b.dense_A.cols() == 2));

}
