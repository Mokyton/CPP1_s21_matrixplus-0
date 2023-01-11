#include <gtest/gtest.h>

#include "S21Matrix.h"

TEST(Test_1, EqMatrix) {
  S21Matrix m1(3, 3);
  S21Matrix m2(3, 3);
  S21Matrix m3(5, 5);
  m1.set_from_0_to_future();
  m2.set_from_0_to_future();
  EXPECT_TRUE(m1.EqMatrix(m2));
  EXPECT_FALSE(m1.EqMatrix(m3));
  m1.set_value(2.0000007, 0, 0);
  m2.set_value(2.0000007, 0, 0);
  EXPECT_TRUE(m1.EqMatrix(m2));
  m2.set_value(2.0000001, 0, 0);
  EXPECT_FALSE(m1.EqMatrix(m2));
}

TEST(Test_2, SumMatrix) {
  S21Matrix m1(4, 4);
  S21Matrix m2(4, 4);
  m1.set_from_0_to_future();
  m2.set_from_0_to_future();
  m1.SumMatrix(m2);
  m2.MulNumber(2);
  EXPECT_TRUE(m1 == m2);
}

TEST(Test_3, SubMatrix) {
  S21Matrix m1(4, 4);
  S21Matrix m2(4, 4);
  m1.set_from_0_to_future();
  m2.set_from_0_to_future();
  m1.SubMatrix(m2);
  S21Matrix m3(4, 4);
  EXPECT_TRUE(m1 == m3);
}

TEST(Test_4, mult_number) {
  S21Matrix m1(2, 4);
  m1.set_from_0_to_future();
  m1.MulNumber(2);
  double tmp[2][4] = {{0, 2, 4, 6}, {8, 10, 12, 14}};
  for (int i = 0; i < m1.get_rows(); i++) {
    for (int j = 0; j < m1.get_cols(); j++) {
      EXPECT_DOUBLE_EQ(m1(i, j), tmp[i][j]);
    }
  }
}

TEST(Test_5, mult_matrix) {
  S21Matrix m1 = S21Matrix(3, 3);
  S21Matrix m2 = S21Matrix(3, 3);
  m1.set_from_0_to_future();
  m2.set_from_0_to_future();
  m1.MulMatrix(m2);
  S21Matrix m3(3, 3);
  m3.set_value(15, 0, 0);
  m3.set_value(18, 0, 1);
  m3.set_value(21, 0, 2);
  m3.set_value(42, 1, 0);
  m3.set_value(54, 1, 1);
  m3.set_value(66, 1, 2);
  m3.set_value(69, 2, 0);
  m3.set_value(90, 2, 1);
  m3.set_value(111, 2, 2);

  for (int i = 0; i < m3.get_rows(); i++) {
    for (int j = 0; j < m1.get_cols(); j++) {
      EXPECT_DOUBLE_EQ(m1(i, j), m3(i, j));
    }
  }
}

TEST(Test_5, transpose) {
  S21Matrix m1(3, 3);
  m1.set_from_0_to_future();
  S21Matrix m2 = m1.Transpose();
  double tmp[3][3] = {{0, 3, 6}, {1, 4, 7}, {2, 5, 8}};
  for (int i = 0; i < m1.get_rows(); i++) {
    for (int j = 0; j < m1.get_cols(); j++) {
      EXPECT_DOUBLE_EQ(m2(i, j), tmp[i][j]);
    }
  }
}

TEST(Test_6, calc_complements) {
  S21Matrix m1(3, 3);
  m1.set_from_0_to_future();
  S21Matrix m2 = m1.CalcComplements();
  double tmp[3][3] = {{-3, 6, -3}, {6, -12, 6}, {-3, 6, -3}};
  for (int i = 0; i < m2.get_rows(); i++) {
    for (int j = 0; j < m2.get_cols(); j++) {
      EXPECT_DOUBLE_EQ(m2(i, j), tmp[i][j]);
    }
  }
}

TEST(Test_7, determinat) {
  S21Matrix m1(3, 3), m2(3, 3);
  m2.set_from_0_to_future();
  double tmp = m2.Determinant();
  EXPECT_DOUBLE_EQ(tmp, 0);
  m1.set_value(98, 0, 0);
  m1.set_value(52, 0, 1);
  m1.set_value(50, 0, 2);
  m1.set_value(30, 1, 0);
  m1.set_value(18, 1, 1);
  m1.set_value(34, 1, 2);
  m1.set_value(1, 2, 0);
  m1.set_value(2, 2, 1);
  m1.set_value(3, 2, 2);
  tmp = m1.Determinant();
  EXPECT_DOUBLE_EQ(tmp, -2184);
}

TEST(Test_8, inverse_matrix) {
  S21Matrix m1(3, 3);
  m1.set_value(98, 0, 0);
  m1.set_value(52, 0, 1);
  m1.set_value(50, 0, 2);
  m1.set_value(30, 1, 0);
  m1.set_value(18, 1, 1);
  m1.set_value(34, 1, 2);
  m1.set_value(1, 2, 0);
  m1.set_value(2, 2, 1);
  m1.set_value(3, 2, 2);
  S21Matrix m2 = m1.InverseMatrix();
  S21Matrix m3(3, 3);
  m3.set_value(0.006410, 0, 0);
  m3.set_value(0.025641, 0, 1);
  m3.set_value(-0.397436, 0, 2);
  m3.set_value(0.025641, 1, 0);
  m3.set_value(-0.111722, 1, 1);
  m3.set_value(0.838828, 1, 2);
  m3.set_value(-0.019231, 2, 0);
  m3.set_value(0.065934, 2, 1);
  m3.set_value(-0.093407, 2, 2);
  for (int i = 0; i < m2.get_rows(); i++) {
    for (int j = 0; j < m2.get_cols(); j++) {
      double val_1 = m2(i, j);
      double val_2 = m3(i, j);
      if (val_1 == val_2) {
        EXPECT_TRUE(1);
      } else {
        EXPECT_FALSE(false);
      }
    }
  }
}

TEST(Test_9, operator_sum) {
  S21Matrix m1(4, 4);
  S21Matrix m2(4, 4);
  S21Matrix m3;
  m1.set_from_0_to_future();
  m2.set_from_0_to_future();
  m3 = m1 + m2;
  m2 *= 2;
  EXPECT_TRUE(m3 == m2);
}

TEST(Test_10, operator_sub) {
  S21Matrix m1(4, 4);
  S21Matrix m2(4, 4);
  S21Matrix m3;
  m1.set_from_0_to_future();
  m2.set_from_0_to_future();
  m3 = m1 - m2;
  m2.zero();
  EXPECT_TRUE(m3 == m2);
}

TEST(Test_11, oper_sub_sum_mult_eq) {
  S21Matrix m1(4, 4);
  S21Matrix m2(4, 4);
  S21Matrix m3(4, 4);
  m1.set_from_0_to_future();
  m2.set_from_0_to_future();
  m3.set_from_0_to_future();
  m1 += m2;
  m3 *= 2;
  EXPECT_TRUE(m1 == m3);
  m3 -= m1;
  m1.zero();
  EXPECT_TRUE(m1 == m3);
}

TEST(Test_12, craze_mega_test) {
  S21Matrix m1(4, 4);
  m1.set_from_0_to_future();
  S21Matrix m2 = std::move(m1);
  S21Matrix m3(m2);
  EXPECT_TRUE(m2 == m3);
}

TEST(Test_13, setters) {
  S21Matrix m1(3, 3);
  m1.set_from_0_to_future();
  double cringe[2][3] = {{0, 1, 2}, {3, 4, 5}};
  m1.set_rows(2);
  for (int i = 0; i < m1.get_rows(); i++) {
    for (int j = 0; j < m1.get_cols(); j++) {
      EXPECT_EQ(m1(i, j), cringe[i][j]);
    }
  }

  double oh_no_cringe[4][3] = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}, {0, 0, 0}};
  S21Matrix m2(3, 3);
  m2.set_from_0_to_future();
  m2.set_rows(4);
  for (int i = 0; i < m2.get_rows(); i++) {
    for (int j = 0; j < m2.get_cols(); j++) {
      EXPECT_EQ(m2(i, j), oh_no_cringe[i][j]);
    }
  }

  S21Matrix m3(2, 2);
  m3.set_from_0_to_future();
  m3.set_cols(1);
  double new_cring[2][1] = {{0}, {2}};
  for (int i = 0; i < m3.get_rows(); i++) {
    for (int j = 0; j < m3.get_cols(); j++) {
      EXPECT_EQ(m3(i, j), new_cring[i][j]);
    }
  }
  S21Matrix m4(2, 2);
  m4.set_from_0_to_future();
  m4.set_cols(3);
  double crg[2][3] = {{0, 1, 0}, {2, 3, 0}};

  for (int i = 0; i < m4.get_rows(); i++) {
    for (int j = 0; j < m4.get_cols(); j++) {
      EXPECT_EQ(m4(i, j), crg[i][j]);
    }
  }
  S21Matrix m5(3, 3);
  m5.set_from_0_to_future();
  double **new_star = m5.get_p();
  EXPECT_EQ(new_star[1][1], 4);
  m5.set_value(0, 1, 1);
  double **blond = m5.get_p();
  EXPECT_EQ(blond[1][1], 0);
}

TEST(Test_14, lucky_me) {
  S21Matrix m1(3, 3);
  S21Matrix m2(3, 3);
  m1.set_from_0_to_future();
  m2.set_from_0_to_future();
  m1 = m2;
  EXPECT_TRUE(m1 == m2);
}

TEST(Test_15, girls_with_green_eyes) {
  S21Matrix m1(3, 3);
  S21Matrix m2(3, 3);
  S21Matrix m3(3, 3);
  m1.set_from_0_to_future();
  m2.set_from_0_to_future();
  S21Matrix m4 = m1 * m2;
  m3.set_value(15, 0, 0);
  m3.set_value(18, 0, 1);
  m3.set_value(21, 0, 2);
  m3.set_value(42, 1, 0);
  m3.set_value(54, 1, 1);
  m3.set_value(66, 1, 2);
  m3.set_value(69, 2, 0);
  m3.set_value(90, 2, 1);
  m3.set_value(111, 2, 2);
  EXPECT_TRUE(m3 == m4);
  m1 *= m2;
  EXPECT_TRUE(m1 == m4);
  S21Matrix m5(m2);
  m5 = m5 * 2;
  m2 += m2;
  EXPECT_TRUE(m2 == m5);
}

// TEST(Test_16, DORA) {
// S21Matrix m1(3,3);
// S21Matrix m2(2, 2);
// EXPECT_THROW(m1 + m2, std::out_of_range);
// EXPECT_THROW(m1 - m2, std::out_of_range);
// EXPECT_THROW(m1 * m2, std::out_of_range);
// S21Matrix m3(3, 2);
// EXPECT_THROW(m3.CalcComplements(), std::out_of_range);
// EXPECT_THROW(m3.Determinant(), std::out_of_range);
// S21Matrix m4(3, 3);
// m4.set_from_0_to_future();
// EXPECT_THROW(m4.InverseMatrix(), std::out_of_range);
// EXPECT_THROW(m4(4, 3), std::out_of_range);
// }

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}