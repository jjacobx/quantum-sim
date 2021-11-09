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
}
