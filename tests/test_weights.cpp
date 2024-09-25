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


TEST_CASE( "Check Weight Properties", "[require]" ){
    //Vaidya, John, DikinLS
    SparseVaidyaWalk vaidya_sparse(0.5);
    SparseDikinLSWalk dikinls_sparse(1.0, 1, 0.001, 0.01, 20000);
    SparseJohnWalk john_sparse(0.5, 1, 1e-5, 10000);

    VectorXd simplex_start (3);
    simplex_start << 0.33, 0.34, 0.33;
    SparseMatrixXd w = dikinls_sparse.generateWeight(simplex_start, simplex.A, simplex.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(2, 0.01));
    w = john_sparse.generateWeight(simplex_start, simplex.A, simplex.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(3, 0.01));
    w = vaidya_sparse.generateWeight(simplex_start, simplex.A, simplex.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(4, 0.01));

    VectorXd hc_start (6);
    hc_start << 0, 0, 1, 1, 1, 1;
    w = dikinls_sparse.generateWeight(hc_start, hc.A, hc.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(2, 0.01));
    w = john_sparse.generateWeight(hc_start, hc.A, hc.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(3, 0.01));
    w = vaidya_sparse.generateWeight(hc_start, hc.A, hc.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(4, 0.01));

    VectorXd birk_start (4);
    birk_start << 0.5, 0.5, 0.5, 0.5;
    w = dikinls_sparse.generateWeight(birk_start, birk.A, birk.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(1, 0.01));
    w = john_sparse.generateWeight(birk_start, birk.A, birk.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(1.5, 0.01));
    w = vaidya_sparse.generateWeight(birk_start, birk.A, birk.k);
    REQUIRE_THAT(w.diagonal().sum(), Catch::Matchers::WithinAbs(2, 0.01));

    FacialReduction fr;
    DenseCenter dc;
    res simplex_dense = fr.reduce(simplex.A, simplex.b, simplex.k, false);
    res hc_dense = fr.reduce(hc.A, hc.b, hc.k, false);
    res birk_dense = fr.reduce(birk.A, birk.b, birk.k, false);
    VectorXd sd_x = dc.getInitialPoint(simplex_dense.dense_A, simplex_dense.dense_b);
    VectorXd hc_x = dc.getInitialPoint(hc_dense.dense_A, hc_dense.dense_b);
    VectorXd birk_x = dc.getInitialPoint(birk_dense.dense_A, birk_dense.dense_b);

    JohnWalk john(0.5, 1, 0.001, 10000);
    DikinLSWalk dikinls(0.5, 1, 0.001, 0.01, 10000);
    VaidyaWalk vaidya(0.5);

    john.generateWeight(sd_x, simplex_dense.dense_A, simplex_dense.dense_b);
    dikinls.generateWeight(sd_x, simplex_dense.dense_A, simplex_dense.dense_b);
    vaidya.generateWeight(sd_x, simplex_dense.dense_A, simplex_dense.dense_b);
    double dw, jw, vw;
    dw = dikinls.weights.diagonal().sum();
    jw = john.weights.diagonal().sum();
    vw = vaidya.weights.diagonal().sum();

    REQUIRE_THAT(dw, Catch::Matchers::WithinAbs(2, 0.01));
    REQUIRE_THAT(jw, Catch::Matchers::WithinAbs(3, 0.01));
    REQUIRE_THAT(vw, Catch::Matchers::WithinAbs(4, 0.01));

    john.generateWeight(hc_x, hc_dense.dense_A, hc_dense.dense_b);
    dikinls.generateWeight(hc_x, hc_dense.dense_A, hc_dense.dense_b);
    vaidya.generateWeight(hc_x, hc_dense.dense_A, hc_dense.dense_b);
    dw = dikinls.weights.diagonal().sum();
    jw = john.weights.diagonal().sum();
    vw = vaidya.weights.diagonal().sum();

    REQUIRE_THAT(dw, Catch::Matchers::WithinAbs(2, 0.01));
    REQUIRE_THAT(jw, Catch::Matchers::WithinAbs(3, 0.01));
    REQUIRE_THAT(vw, Catch::Matchers::WithinAbs(4, 0.01));

    john.generateWeight(birk_x, birk_dense.dense_A, birk_dense.dense_b);
    dikinls.generateWeight(birk_x, birk_dense.dense_A, birk_dense.dense_b);
    vaidya.generateWeight(birk_x, birk_dense.dense_A, birk_dense.dense_b);
    dw = dikinls.weights.diagonal().sum();
    jw = john.weights.diagonal().sum();
    vw = vaidya.weights.diagonal().sum();

    REQUIRE_THAT(dw, Catch::Matchers::WithinAbs(1, 0.01));
    REQUIRE_THAT(jw, Catch::Matchers::WithinAbs(1.5, 0.01));
    REQUIRE_THAT(vw, Catch::Matchers::WithinAbs(2, 0.01));

}