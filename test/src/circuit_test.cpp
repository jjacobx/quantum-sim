#include <limits>
#include "gtest/gtest.h"
#include "circuit.h"

class GateTest : public::testing::Test {
  protected:
    virtual void SetUp() {
      // code before each test
    }

    virtual void TearDown() {
      // code after each test
    }
};

TEST_F(GateTest, gate_comparison) {
  EXPECT_EQ(Ops::I, Ops::I);
  EXPECT_NE(Ops::X, Ops::Y);
}

TEST_F(GateTest, gate_comparison_invalid) {
  EXPECT_THROW(Ops::I == Ops::CNOT, invalid_argument);
  EXPECT_THROW(Ops::CNOT == Ops::CCNOT, invalid_argument);
}

TEST_F(GateTest, gate_creation_sparse_dense) {
  DMatrix m_dense {{1, 0}, {0, 1}};
  SMatrix m_sparse = m_dense.sparseView();
  EXPECT_EQ(Gate(m_dense), Gate(m_sparse));
}

TEST_F(GateTest, gate_creation_invalid_dimensions) {
  DMatrix m_invalid_1 {{1, 2, 3}, {4, 5, 6}};
  DMatrix m_invalid_2 {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  EXPECT_THROW(Gate g_invalid_1 = Gate(m_invalid_1), invalid_argument);
  EXPECT_THROW(Gate g_invalid_2 = Gate(m_invalid_2), invalid_argument);
}

TEST_F(GateTest, gate_creation_non_unitary) {
  DMatrix m_non_unitary {{1, 2}, {3, 4}};
  EXPECT_THROW(Gate g_non_unitary = Gate(m_non_unitary), invalid_argument);
}

TEST_F(GateTest, gate_multiplication_dense) {
  EXPECT_EQ(Ops::X * Ops::X, Ops::I);
  EXPECT_EQ(Ops::Y * Ops::Y, Ops::I);
  EXPECT_EQ(Ops::Z * Ops::Z, Ops::I);
  EXPECT_EQ(Ops::H * Ops::H, Ops::I);
  EXPECT_EQ(Ops::H * Ops::H * Ops::H, Ops::H);
}

TEST_F(GateTest, gate_multiplication_sparse) {
  EXPECT_EQ(Ops::CNOT * Ops::CNOT, Ops::Id(2));
  EXPECT_EQ(Ops::CCNOT * Ops::CCNOT, Ops::Id(3));
}

TEST_F(GateTest, gate_multiplication_invalid) {
  EXPECT_THROW(Ops::I * Ops::CNOT, invalid_argument);
  EXPECT_THROW(Ops::CNOT * Ops::CCNOT, invalid_argument);
}

TEST_F(GateTest, gate_kronecker_dense) {
  DMatrix m_kronecker {
    {0.5,  0.5,  0.5,  0.5}, 
    {0.5, -0.5,  0.5, -0.5}, 
    {0.5,  0.5, -0.5, -0.5}, 
    {0.5, -0.5, -0.5,  0.5}
  };
  EXPECT_EQ(Ops::H & Ops::H, Gate(m_kronecker));
  EXPECT_EQ(Ops::I & Ops::I & Ops::I, Ops::Id(3));
}

TEST_F(GateTest, gate_kronecker_sparse) {
  EXPECT_EQ(Ops::Id(2) & Ops::Id(2), Ops::Id(4));
  EXPECT_EQ(Ops::Id(2) & Ops::Id(3) & Ops::Id(2), Ops::Id(7));
}

TEST_F(GateTest, gate_adjoint) {
  DMatrix m_anti_hermitian {{0, 1}, {-1, 0}};
  EXPECT_EQ(Ops::X.adjoint(), Ops::X);
  EXPECT_EQ(Gate(m_anti_hermitian).adjoint(), -Gate(m_anti_hermitian));
  EXPECT_EQ(Ops::CCNOT.adjoint(), Ops::CCNOT);
}
