#include "circuit.h"
#include <unsupported/Eigen/KroneckerProduct>

ostream& operator<<(ostream& os, Gate const& g) {
  return g.print(os);
}

MatrixRep Gate::rep() const {
  if (holds_alternative<MatrixXcd>(this->op))
    return MatrixRep::dense;
  else return MatrixRep::sparse;
}

ostream& Gate::print(ostream& os) const {
  if (this->rep() == MatrixRep::dense)
    return os << get<DMatrix>(this->get_op());
  else return os << get<SMatrix>(this->get_op());
}

bool Gate::operator==(const Gate& rhs) const {
  if (this->get_dim() != rhs.get_dim())
    throw invalid_argument("Cannot compare gates with different dimensions");

  return visit([&](auto&& m1) -> bool { 
    return visit([&](auto&& m2) -> bool {
      return (m1 - m2).norm() < 0.01;
    }, rhs.get_op());
  }, this->get_op());
}

bool Gate::operator!=(const Gate& rhs) const {
  return !(*this == rhs);
}

Gate Gate::operator*(const Gate& rhs) const {
  if (this->get_dim() != rhs.get_dim())
    throw invalid_argument("Cannot multiply gates with different dimensions");

  MatrixVar res = visit([&](auto&& m1) -> MatrixVar { 
    return visit([&](auto&& m2) -> MatrixVar {
      return (m1 * m2).eval();
    }, rhs.get_op());
  }, this->get_op());
  return Gate(res);
}

Gate Gate::operator&(const Gate& rhs) const {
  MatrixVar res;
  if (this->rep() == MatrixRep::dense && rhs.rep() == MatrixRep::dense) {
    res = (DMatrix)KroneckerProduct<DMatrix, DMatrix>(get<DMatrix>(this->get_op()), get<DMatrix>(rhs.get_op()));
  } else if (this->rep() == MatrixRep::dense && rhs.rep() == MatrixRep::sparse) {
    res = (SMatrix)KroneckerProductSparse<DMatrix, SMatrix>(get<DMatrix>(this->get_op()), get<SMatrix>(rhs.get_op()));
  } else if (this->rep() == MatrixRep::sparse && rhs.rep() == MatrixRep::dense) {
    res = (SMatrix)KroneckerProductSparse<SMatrix, DMatrix>(get<SMatrix>(this->get_op()), get<DMatrix>(rhs.get_op()));
  } else {
    res = (SMatrix)KroneckerProductSparse<SMatrix, SMatrix>(get<SMatrix>(this->get_op()), get<SMatrix>(rhs.get_op()));
  }

  return Gate(res);
}

Gate Gate::operator^(int n) const {
  if (n < 0)
    throw invalid_argument("n needs to be a non-negative integer");
  if (n == 0)
    return Ops::Id(this->qubits.size());
  else
    return *this * (*this ^ (n - 1));
}

Gate Gate::operator-() const {
  MatrixVar res = visit([&](auto&& m) -> MatrixVar { 
    return (-m).eval();
  }, this->get_op());

  return Gate(res);
}

Gate Gate::adjoint() const {
  MatrixVar res = visit([&](auto&& m) -> MatrixVar { 
    return (decltype(m))m.adjoint().eval();
  }, this->get_op());
  
  return Gate(res);
}

void Gate::assign(MatrixVar op_) {
  if (this->rep() == MatrixRep::dense)
    assign(std::get<DMatrix>(op_));
  else  assign(std::get<SMatrix>(op_));
}

void Gate::check_unitarity() {
  bool is_unitary = visit([&](auto&& m) -> bool { 
    SMatrix id(this->dim, this->dim);
    id.setIdentity();
    return (m * m.adjoint() - id).norm() < 0.01;
  }, this->get_op());

  if (!is_unitary)
    throw invalid_argument("The operator needs to be a unitary matrix");
}
