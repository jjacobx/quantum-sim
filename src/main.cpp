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

  SMatrixcd res = KroneckerProductSparse<SMatrixcd, Matrix2cd>(Ops::CNOT, Ops::H);
  cout << "CNOT o H = " << endl << res << endl;
}
