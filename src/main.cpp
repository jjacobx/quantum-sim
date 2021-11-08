#include "circuit.h"
#include <iostream>
#include <unsupported/Eigen/KroneckerProduct>


using namespace std;

int main()
{
  cout << "I = " << endl << Ops::I << endl;
  cout << "X = " << endl << Ops::X << endl;
  cout << "Y = " << endl << Ops::Y << endl;
  cout << "Z = " << endl << Ops::Z << endl;
  cout << "H = " << endl << Ops::H << endl;
  cout << "S = " << endl << Ops::S << endl;
  cout << "T = " << endl << Ops::T << endl;
  cout << "CNOT = " << endl << Ops::CNOT << endl;
  cout << "CCNOT = " << endl << Ops::CCNOT << endl;

  cout << "Rx(PI/4) = " << endl << Ops::Rx(M_PI / 4.0) << endl;
  cout << "Ry(PI/4) = " << endl << Ops::Ry(M_PI / 4.0) << endl;
  cout << "Rz(PI/4) = " << endl << Ops::Rz(M_PI / 4.0) << endl;

  SMatrixcd res = KroneckerProductSparse<SMatrixcd, Matrix2cd>(Ops::CNOT, Ops::H);
  cout << "CNOT o H = " << endl << res << endl;
  cout << (int)sqrt(Ops::X.size()) << endl;
  cout << Ops::Expi(Ops::X, M_PI / 2.0) << endl;
}
