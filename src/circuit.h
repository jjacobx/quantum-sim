#pragma once
#include <complex>
#include <iostream>
#include <stdexcept>
#include <variant>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;
using namespace std::complex_literals;

typedef MatrixXcd DMatrix;
typedef SparseMatrix<dcomplex> SMatrix;
typedef variant<DMatrix, SMatrix> MatrixVar;

enum class MatrixRep { dense, sparse };

class Gate {
  public: 
    Gate(MatrixVar op_) : op(op_) { assign(op_); }

    MatrixVar get_op() const { return this->op; }
    int get_dim() const { return this->dim; }

    MatrixRep rep() const;                  // representation type
    ostream& print(ostream& os) const;      // printing
    Gate operator*(const Gate& rhs) const;  // multiplication
    Gate operator&(const Gate& rhs) const;  // tensor product
    Gate operator^(int n) const;            // integer powers

  private:
    MatrixVar op;
    int dim;
    vector<int> qubits;

    template<class Derived>
    void assign(EigenBase<Derived>& op_) {
      if (op_.rows() != op_.cols())
        throw std::invalid_argument("Operator needs to be a square matrix");
      dim = op_.rows();
      if (dim == 1 || (dim & (dim - 1)))
        throw std::invalid_argument("Dimensions of the operator need to be positive integer powers of 2");
      qubits.assign((int)log2(dim), -1);
    }

    void assign(MatrixVar op_);
};

namespace Ops {
  inline Gate I {
    Matrix2cd::Identity()
  };
  inline Gate X {
    DMatrix({
      {0.0, 1.0},
      {1.0, 0.0}
    })
  };
  inline Gate Y { 
    DMatrix({
     {0.0,  -1.0i},
     {1.0i,   0.0}
    })
  };
  inline Gate Z { 
    DMatrix({
     {1.0,   0.0},
     {0.0i, -1.0}
    })
  };
  inline Gate H { 
    DMatrix({
     {1 / sqrt(2),  1 / sqrt(2)},
     {1 / sqrt(2), -1 / sqrt(2)}
    })
  };
  inline Gate S { 
    DMatrix({
     {1.0,  0.0},
     {0.0, 1.0i}
    })
  };
  inline Gate T { 
    DMatrix({
     {1.0,                    0.0},
     {0.0, exp(1.0i * M_PI / 4.0)}
    })
  };
  inline Gate CNOT { 
    (SMatrix)DMatrix({
      {1.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0}, 
      {0.0, 0.0, 0.0, 1.0}, 
      {0.0, 0.0, 1.0, 0.0}
    }).sparseView() 
  };
  inline Gate CCNOT { 
    (SMatrix)DMatrix({
     {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
     {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
     {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
     {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, 
     {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, 
     {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0}, 
     {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}, 
     {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0}
    }).sparseView() 
  };

  inline Gate Id(int n) {
    if (n == 2) return I;
    SMatrix id(n, n);
    id.setIdentity();
    return Gate(id);
  }

  inline Gate Rx(double theta) {
    return(Gate(DMatrix({
      {  1.0 * cos(theta / 2.0), -1.0i * sin(theta / 2.0)},
      {-1.0i * sin(theta / 2.0),   1.0 * cos(theta / 2.0)}
    })));
  }

  inline Gate Ry(double theta) {
    return(Gate(DMatrix({
      {cos(theta / 2.0), -sin(theta / 2.0)},
      {sin(theta / 2.0),  cos(theta / 2.0)}
    })));
  }

  inline Gate Rz(double theta) {
    return(Gate(DMatrix({
      {exp(-1.0i * theta / 2.0),                      0.0},
      {                     0.0, exp(-1.0i * theta / 2.0)}
    })));
  }
}

namespace qol {
  inline ostream& operator<<(ostream& os, Gate g) {
    return g.print(os);
  }
}
