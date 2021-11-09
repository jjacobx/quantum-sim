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

TEST_F(GateTest, gate_multiplication) {
  EXPECT_EQ(Ops::X * Ops::X, Ops::I);
  EXPECT_EQ(Ops::Y * Ops::Y, Ops::I);
  EXPECT_EQ(Ops::Z * Ops::Z, Ops::I);
  EXPECT_EQ(Ops::H * Ops::H, Ops::I);
}
