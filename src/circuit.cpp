#include "circuit.h"

using namespace std::complex_literals;

Matrix2cd Ops::I = Matrix2cd::Identity();
Matrix2cd Ops::X {
     {0.0, 1.0},
     {1.0, 0.0}
};
Matrix2cd Ops::Y {
     {0.0,  -1.0i},
     {1.0i,   0.0}
};
Matrix2cd Ops::Z {
     {1.0,   0.0},
     {0.0i, -1.0}
};
Matrix2cd Ops::H {
     {1 / sqrt(2),  1 / sqrt(2)},
     {1 / sqrt(2), -1 / sqrt(2)}
};
Matrix2cd Ops::S {
     {1.0,  0.0},
     {0.0, 1.0i}
};
Matrix2cd Ops::T {
     {1.0,                    0.0},
     {0.0, exp(1.0i * M_PI / 4.0)}
};

SMatrixcd cnot() {
    SparseMatrix<dcomplex> cnot(4, 4);
    cnot.coeffRef(0, 0) = 1;
    cnot.coeffRef(1, 1) = 1;
    cnot.coeffRef(2, 3) = 1;
    cnot.coeffRef(3, 2) = 1;
    return(cnot);
}
SMatrixcd Ops::CNOT = cnot();

SMatrixcd ccnot() {
    SparseMatrix<dcomplex> ccnot(8, 8);
    for (int i = 0; i < 6; i++)
        ccnot.coeffRef(i, i) = 1;
    ccnot.coeffRef(6, 7) = 1;
    ccnot.coeffRef(7, 6) = 1;
    return(ccnot);
}
SMatrixcd Ops::CCNOT = ccnot();
