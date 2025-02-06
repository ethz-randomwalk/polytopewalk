#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "utils/FullWalkRun.hpp"
#include <cstring>


double uniform_dist(const VectorXd& x) { return 1.0; };

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

TEST_CASE( "Test All Dense Combinations", "[require]" ){
    JohnWalk john(0.5, 0, uniform_dist, 1, 1e-5, 1000);
    DikinLSWalk dikinls(3.0, 0, uniform_dist, 1, 0.001, 0.01, 100);
    VaidyaWalk vaidya(0.5, 0, uniform_dist);
    DikinWalk dikin(0.5, 0, uniform_dist);
    BallWalk ball(0.5);
    HitAndRun hitrun(0.5, 1, 0.001);
    DenseCenter dc;
    FacialReduction fr; 

    MatrixXd walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &john, &fr, &dc);
    walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &dikinls, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &vaidya, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &dikin, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &ball, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(simplex.A, simplex.b, simplex.k, 100, &hitrun, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);

    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &john, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &dikinls, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &vaidya, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &dikin, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &ball, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(hc.A, hc.b, hc.k, 100, &hitrun, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);

    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &john, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &dikinls, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &vaidya, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &dikin, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &ball, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
    walk_res = denseFullWalkRun(birk.A, birk.b, birk.k, 100, &hitrun, &fr, &dc);
    REQUIRE(walk_res.rows() == 100);
}