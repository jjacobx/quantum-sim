#include "circuit.h"
#include <unsupported/Eigen/KroneckerProduct>

MatrixRep Gate::rep() const {
  if (std::holds_alternative<MatrixXcd>(this->op))
    return MatrixRep::dense;
  else return MatrixRep::sparse;
}

ostream& Gate::print(ostream& os) const {
  if (this->rep() == MatrixRep::dense)
    return os << get<MatrixXcd>(this->get_op());
  else return os << get<SMatrix>(this->get_op());
}

Gate Gate::operator*(const Gate& rhs) const {
  MatrixVar res = std::visit([&](auto&& m1) -> MatrixVar { 
    return std::visit([&](auto&& m2) -> MatrixVar {
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
    return Ops::Id(this->dim);
  else
    return *this * (*this ^ (n - 1));
}

void Gate::assign(MatrixVar op_) {
  if (this->rep() == MatrixRep::dense)
    assign(std::get<DMatrix>(op_));
  else  assign(std::get<SMatrix>(op_));
}
