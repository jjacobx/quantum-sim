#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#pragma once

using namespace Eigen;
typedef SparseMatrix<dcomplex> SMatrixcd;

class Ops {
  public: 
    static Matrix2cd I, X, Y, Z, H, S, T;
    static SMatrixcd CNOT, CCNOT;
};
