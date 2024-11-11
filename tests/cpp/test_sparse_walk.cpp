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
    SparseJohnWalk john(0.5, 2);
    SparseDikinLSWalk dikinls(3.0, 2);
    SparseVaidyaWalk vaidya(0.5, 2);
    SparseDikinWalk dikin(0.5, 2);
    SparseBallWalk ball(0.5, 2);
    SparseHitAndRun hitrun(0.5, 2);
    SparseCenter sc;
    FacialReduction fr; 

    MatrixXd walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &john, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &dikinls, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &vaidya, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &dikin, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &ball, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &hitrun, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);

    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &john, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &dikinls, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &vaidya, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &dikin, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &ball, &fr, &sc, 1);;
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(hc.A, hc.b, hc.k, 100, &hitrun, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);

    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &john, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &dikinls, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &vaidya, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &dikin, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &ball, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
    walk_res = sparseFullWalkRun(birk.A, birk.b, birk.k, 100, &hitrun, &fr, &sc, 1);
    REQUIRE(walk_res.rows() == 100);
}