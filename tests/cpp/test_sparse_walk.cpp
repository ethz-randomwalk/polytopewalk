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

TEST_CASE( "Test All Sparse Combinations", "[require]" ){
    SparseJohnWalk john(0.5);
    SparseDikinLSWalk dikinls(3.0);
    SparseVaidyaWalk vaidya(0.5);
    SparseDikinWalk dikin(0.5);
    SparseBallWalk ball(0.5);
    SparseHitAndRun hitrun(0.5);
    SparseCenter sc;
    FacialReduction fr; 

    MatrixXd walk_res = sparseFullWalkRun(100, simplex.A, simplex.b, simplex.k, &john, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, simplex.A, simplex.b, simplex.k, &dikinls, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, simplex.A, simplex.b, simplex.k, &vaidya, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, simplex.A, simplex.b, simplex.k, &dikin, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, simplex.A, simplex.b, simplex.k, &ball, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, simplex.A, simplex.b, simplex.k, &hitrun, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);

    walk_res = sparseFullWalkRun(100, hc.A, hc.b, hc.k, &john, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, hc.A, hc.b, hc.k, &dikinls, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, hc.A, hc.b, hc.k, &vaidya, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, hc.A, hc.b, hc.k, &dikin, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, hc.A, hc.b, hc.k, &ball, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, hc.A, hc.b, hc.k, &hitrun, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);

    walk_res = sparseFullWalkRun(100, birk.A, birk.b, birk.k, &john, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, birk.A, birk.b, birk.k, &dikinls, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, birk.A, birk.b, birk.k, &vaidya, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, birk.A, birk.b, birk.k, &dikin, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, birk.A, birk.b, birk.k, &ball, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
    walk_res = sparseFullWalkRun(100, birk.A, birk.b, birk.k, &hitrun, &fr, &sc, 10, 10, 1220);
    REQUIRE(walk_res.rows() == 9);
}