#ifndef SRC_S21MATRIX_H_
#define SRC_S21MATRIX_H_

#include <math.h>

#include <cstring>
#include <iostream>
//#include <stdio.h>
class S21Matrix {
 private:
  // attributes
  int _rows, _cols;  // rows and columns attributes
  double** _p;
  // additional methods
  void copy_obj(const S21Matrix& o);
  void copy_matrix(const S21Matrix& o, int i);
  void create_matrix();
  void remove_matrix();

 public:
  S21Matrix();                    // default constructor
  S21Matrix(int rows, int cols);  // parameterized constructor
  S21Matrix(const S21Matrix& o);  // copy cnstructor
  S21Matrix(S21Matrix&& o);       // move cnstructor
  ~S21Matrix();                   // destructor

  // some operators overloads
  S21Matrix& operator=(const S21Matrix& o);  // assignment operator overload
  double& operator()(int row, int col);      // index operator overload
  S21Matrix& operator+=(const S21Matrix& o);
  S21Matrix operator+(const S21Matrix& o);
  S21Matrix& operator-=(const S21Matrix& o);
  S21Matrix operator-(const S21Matrix& o);
  S21Matrix operator*(const double num);
  S21Matrix& operator*=(const double num);
  S21Matrix operator*(const S21Matrix& o);
  S21Matrix& operator*=(const S21Matrix& o);
  bool operator==(const S21Matrix& o);

  // some public methods
  bool EqMatrix(const S21Matrix& o);
  void SumMatrix(const S21Matrix& o);
  void SubMatrix(const S21Matrix& o);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& o);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  S21Matrix looking_for_minor(int rows, int columns);
  double Determinant();
  S21Matrix InverseMatrix();

  // some getters
  int get_rows();
  int get_cols();
  double** get_p();

  // some setters
  void set_rows(int rows);
  void set_cols(int cols);
  void set_value(double value, int rows, int cols);
  void set_from_0_to_future();

  // additional methods
  void m_out();
  void zero();
};

#endif  // SRC_S21MATRIX_H_