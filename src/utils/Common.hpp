#ifndef COMMON1_HPP
#define COMMON1_HPP

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <unsupported/Eigen/MatrixFunctions>
typedef Eigen::SparseMatrix<double> SparseMatrixXd;
typedef Eigen::Triplet<double> T;

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <cstring>
// use extern "C" to the C++ compiler 
// to treat the code inside as if it were C code
#ifdef __cplusplus
extern "C" {
#endif

#include <glpk.h>

#ifdef __cplusplus
}
#endif

using namespace Eigen;
using namespace std;

#endif