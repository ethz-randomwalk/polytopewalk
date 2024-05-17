#ifndef COMMON_HPP
#define COMMON_HPP

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

using namespace Eigen;
using namespace std;

#endif