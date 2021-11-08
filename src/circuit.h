#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#pragma once

using namespace Eigen;
using namespace std::complex_literals;

typedef SparseMatrix<dcomplex> SMatrixcd;
typedef EigenBase<dcomplex> EigenBasecd;

namespace Ops {
  inline Matrix2cd I = Matrix2cd::Identity();
  inline Matrix2cd X {
     {0.0, 1.0},
     {1.0, 0.0}
  };
  inline Matrix2cd Y {
     {0.0,  -1.0i},
     {1.0i,   0.0}
  };
  inline Matrix2cd Z {
     {1.0,   0.0},
     {0.0i, -1.0}
  };
  inline Matrix2cd H {
     {1 / sqrt(2),  1 / sqrt(2)},
     {1 / sqrt(2), -1 / sqrt(2)}
  };
  inline Matrix2cd S {
     {1.0,  0.0},
     {0.0, 1.0i}
  };
  inline Matrix2cd T {
     {1.0,                    0.0},
     {0.0, exp(1.0i * M_PI / 4.0)}
  };
  inline SMatrixcd CNOT = MatrixXcd({
     {1.0, 0.0, 0.0, 0.0},
     {0.0, 1.0, 0.0, 0.0}, 
     {0.0, 0.0, 0.0, 1.0}, 
     {0.0, 0.0, 1.0, 0.0}
  }).sparseView();
  inline SMatrixcd CCNOT = MatrixXcd({
     {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
     {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
     {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
     {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, 
     {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, 
     {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0}, 
     {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}, 
     {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0}, 
  }).sparseView();

  inline Matrix2cd Rx(double theta) {
    return(Matrix2cd({
      {  1.0 * cos(theta / 2.0), -1.0i * sin(theta / 2.0)},
      {-1.0i * sin(theta / 2.0),   1.0 * cos(theta / 2.0)}
    }));
  }

  inline Matrix2cd Ry(double theta) {
    return(Matrix2cd({
      {cos(theta / 2.0), -sin(theta / 2.0)},
      {sin(theta / 2.0),  cos(theta / 2.0)}
    }));
  }

  inline Matrix2cd Rz(double theta) {
    return(Matrix2cd({
      {exp(-1.0i * theta / 2.0),                      0.0},
      {                     0.0, exp(-1.0i * theta / 2.0)}
    }));
  }

  template<class Derived>
  inline typename Derived::PlainObject Expi(EigenBase<Derived>& op_, dcomplex z) {
    int n = op_.rows();
    SMatrixcd id = SMatrixcd(n, n);
    id.setIdentity();

    Derived const& op = op_.derived();
    return(cos(z) * id + 1.0i * sin(z) * id * op);
  }
}
