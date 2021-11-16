#include "circuit.h"
#include <iostream>

using namespace qol;

int main()
{
  cout << "X = " << endl << Ops::X << endl;
  cout << "X * X = " << endl << Ops::X * Ops::X << endl;
  cout << "X ^ 0 = " << endl << (Ops::X ^ 0) << endl;
  cout << "CNOT = " << endl << Ops::CNOT << endl;
  cout << "CNOT ^ 3 = " << endl << (Ops::CNOT ^ 3) << endl;
  cout << "X & CNOT = " << endl << (Ops::X & Ops::CNOT) << endl;
  cout << "X ^ 2 == ID" << endl << ((Ops::X ^ 2) == Ops::I) << endl;
  cout << "Y.adj = " << endl << Ops::Y.adjoint() << endl;
  cout << "Y * Y.adj = " << endl << Ops::Y * Ops::Y.adjoint() << endl;
}
