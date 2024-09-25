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

TEST_CASE( "Check Initialization Algorithm", "[require]" ){
    SparseCenter sc;
    VectorXd simplex_x = sc.getInitialPoint(simplex.A, simplex.b, simplex.k);
    REQUIRE(simplex_x.rows() == 3);
    REQUIRE(simplex_x(0) == Catch::Approx(0.3333333).epsilon(0.01));
    REQUIRE(simplex_x(1) == Catch::Approx(0.3333333).epsilon(0.01));
    REQUIRE(simplex_x(2) == Catch::Approx(0.3333333).epsilon(0.01));

    VectorXd hc_x = sc.getInitialPoint(hc.A, hc.b, hc.k);
    REQUIRE(hc_x.rows() == 6);
    REQUIRE_THAT(hc_x(0), Catch::Matchers::WithinAbs(0, 0.0001));
    REQUIRE_THAT(hc_x(1), Catch::Matchers::WithinAbs(0, 0.0001));
    REQUIRE_THAT(hc_x(2), Catch::Matchers::WithinAbs(1, 0.0001));
    REQUIRE_THAT(hc_x(3), Catch::Matchers::WithinAbs(1, 0.0001));
    REQUIRE_THAT(hc_x(4), Catch::Matchers::WithinAbs(1, 0.0001));
    REQUIRE_THAT(hc_x(5), Catch::Matchers::WithinAbs(1, 0.0001));

    VectorXd birk_x = sc.getInitialPoint(birk.A, birk.b, birk.k);
    REQUIRE(birk_x.rows() == 4);
    REQUIRE_THAT(birk_x(0), Catch::Matchers::WithinAbs(0.5, 0.0001));
    REQUIRE_THAT(birk_x(1), Catch::Matchers::WithinAbs(0.5, 0.0001));
    REQUIRE_THAT(birk_x(2), Catch::Matchers::WithinAbs(0.5, 0.0001));
    REQUIRE_THAT(birk_x(3), Catch::Matchers::WithinAbs(0.5, 0.0001));

    DenseCenter dc; 

    MatrixXd A1 (4, 2);
    A1 << 1, 0, 0, 1, -1, 0, 0, -1;
    VectorXd b1 (4);
    b1 << 1, 1, 0, 0;

    VectorXd center1 = dc.getInitialPoint(A1, b1);
    REQUIRE_THAT(center1(0), Catch::Matchers::WithinAbs(0.5, 0.0001));
    REQUIRE_THAT(center1(1), Catch::Matchers::WithinAbs(0.5, 0.0001));

    MatrixXd A2 (3, 2);
    A2 << -1, 0, 0, -1, 1, 1;
    VectorXd b2 (3);
    b2 << 0, 0, 1;

    VectorXd center2 = dc.getInitialPoint(A2, b2);
    REQUIRE_THAT(center2(0), Catch::Matchers::WithinAbs(0.33, 0.01));
    REQUIRE_THAT(center2(1), Catch::Matchers::WithinAbs(0.33, 0.01));

}