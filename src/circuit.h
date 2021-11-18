#pragma once
#include <complex>
#include <iostream>
#include <stdexcept>
#include <variant>
#include <vector>

/**
 * @file circuit.h
 * @brief Quantum circuit components. 
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;
using namespace std::complex_literals;

typedef MatrixXcd DMatrix;
typedef SparseMatrix<dcomplex> SMatrix;
typedef variant<DMatrix, SMatrix> MatrixVar;

class QMatrix;
ostream& operator<<(ostream& os, QMatrix const& q);

enum class EigenRep { dense, sparse };

class QMatrix {
  public:
    /** @returns The underlying matrix operator as a variant */
    MatrixVar get_op() const { return this->op; }

    EigenRep rep() const;                      /** representation type */
    ostream& print(ostream& os) const;         /** printing */
    bool operator==(const QMatrix& rhs) const; /** comparison */
    bool operator!=(const QMatrix& rhs) const; /** negative comparison */

    // virtual QMatrix operator*(const QMatrix& rhs) const;  /** multiplication */
    // virtual QMatrix operator&(const QMatrix& rhs) const;  /** tensor product */

    friend ostream& operator<<(ostream& os, QMatrix const& q);

  protected:
    MatrixVar op;
};

/**
 * @class Gate
 * @brief A quantum gate represented by either dense or sparse matrix. 
 */
class Gate : public QMatrix {
  public: 
    Gate(MatrixVar op_) { assign(op_); }

    /** @returns The number of columns/rows of a the operator */
    int get_dim() const { return this->dim; }

    Gate operator*(const Gate& rhs) const;  /** multiplication */
    Gate operator&(const Gate& rhs) const;  /** tensor product */
    Gate operator^(int n) const;            /** integer powers */
    Gate operator-() const;                 /** unary negation */
    Gate adjoint() const;

  private:
    int dim;
    vector<int> qubits;

    /**
     * @brief assign class members
     * @param op_ a unitary dense or sparse matrix
     */
    void assign(MatrixVar op_);
    void check_unitarity();

    /** Helper method to match either dense or sparse matrix*/
    template<class Derived>
    void assign(EigenBase<Derived>& op_) {
      if (op_.rows() != op_.cols())
        throw invalid_argument("Operator needs to be a square matrix");
      dim = op_.rows();
      if (dim == 1 || (dim & (dim - 1)))
        throw invalid_argument("Dimensions of the operator need to be positive integer powers of 2");
      qubits.assign((int)log2(dim), -1);

      check_unitarity();
    }

    
};

namespace Ops {
  /** 1 qubit identity gate */
  inline Gate I {
    Matrix2cd::Identity()
  };
  /** 1 qubit Pauli X gate */
  inline Gate X {
    DMatrix({
      {0.0, 1.0},
      {1.0, 0.0}
    })
  };
  /** 1 qubit Pauli Y gate */
  inline Gate Y { 
    DMatrix({
     {0.0,  -1.0i},
     {1.0i,   0.0}
    })
  };
  /** 1 qubit Pauli Z gate */
  inline Gate Z { 
    DMatrix({
     {1.0,   0.0},
     {0.0i, -1.0}
    })
  };
  /** 1 qubit Hadamard gate */
  inline Gate H { 
    DMatrix({
     {1 / sqrt(2),  1 / sqrt(2)},
     {1 / sqrt(2), -1 / sqrt(2)}
    })
  };
  /** 1 qubit phase gate */
  inline Gate S { 
    DMatrix({
     {1.0,  0.0},
     {0.0, 1.0i}
    })
  };
  /** 1 qubit PI/8 gate */
  inline Gate T { 
    DMatrix({
     {1.0,                    0.0},
     {0.0, exp(1.0i * M_PI / 4.0)}
    })
  };
  /** 2 qubit controlled-NOT gate */
  inline Gate CNOT { 
    (SMatrix)DMatrix({
      {1.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0}, 
      {0.0, 0.0, 0.0, 1.0}, 
      {0.0, 0.0, 1.0, 0.0}
    }).sparseView() 
  };
  /** 3 qubit Toffoli gate */
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

  /**
   * @param n_qubits number of qubits to be acted on
   * @return n-qubit identity gate
   */
  inline Gate Id(int n_qubits) {
    int n = 1 << n_qubits;
    if (n == 2) return I;
    SMatrix id(n, n);
    id.setIdentity();
    return Gate(id);
  }

  /**
   * @param theta angle of rotation
   * @return 1-qubit x-axis rotation matrix
   */
  inline Gate Rx(double theta) {
    return(Gate(DMatrix({
      {  1.0 * cos(theta / 2.0), -1.0i * sin(theta / 2.0)},
      {-1.0i * sin(theta / 2.0),   1.0 * cos(theta / 2.0)}
    })));
  }

  /**
   * @param theta angle of rotation
   * @return 1-qubit y-axis rotation matrix
   */
  inline Gate Ry(double theta) {
    return(Gate(DMatrix({
      {cos(theta / 2.0), -sin(theta / 2.0)},
      {sin(theta / 2.0),  cos(theta / 2.0)}
    })));
  }

  /**
   * @param theta angle of rotation
   * @return 1-qubit z-axis rotation matrix
   */
  inline Gate Rz(double theta) {
    return(Gate(DMatrix({
      {exp(-1.0i * theta / 2.0),                      0.0},
      {                     0.0, exp(-1.0i * theta / 2.0)}
    })));
  }
}
