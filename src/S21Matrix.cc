#include "S21Matrix.h"

S21Matrix::S21Matrix() {
  _rows = 3;
  _cols = 3;
  this->create_matrix();
}

S21Matrix::S21Matrix(int rows, int cols) : _rows(rows), _cols(cols) {
  this->create_matrix();
}

S21Matrix::S21Matrix(const S21Matrix& o) : _rows(o._rows), _cols(o._cols) {
  this->create_matrix();
  this->copy_obj(o);
}

S21Matrix::S21Matrix(S21Matrix&& o) {
  if (this != &o) {
    this->_p = o._p;
    this->_rows = o._rows;
    this->_cols = o._cols;
    o._p = NULL;
    o._rows = 0;
    o._cols = 0;
  }
}

S21Matrix::~S21Matrix() {
  this->remove_matrix();
  this->_rows = 0;
  this->_cols = 0;
}

bool S21Matrix::EqMatrix(const S21Matrix& o) {
  bool res = true;
  if (this->_rows != o._rows || this->_cols != o._cols) {
    res = false;
  } else {
    for (auto i = 0; i < _rows; i++) {
      for (auto j = 0; j < _cols; j++) {
        if (fabs(_p[i][j] - o._p[i][j]) > 1e-7 * fabs(o._p[i][j])) {
          res = false;
          break;
        }
        if (!res) break;
      }
    }
  }
  return res;
}

void S21Matrix::SumMatrix(const S21Matrix& o) {
  if (_rows != o._rows || _cols != o._cols) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }
  for (auto i = 0; i < _rows; i++) {
    for (auto j = 0; j < _cols; j++) {
      _p[i][j] += o._p[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& o) {
  if (_rows != o._rows || _cols != o._cols) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }
  for (auto i = 0; i < _rows; i++) {
    for (auto j = 0; j < _cols; j++) {
      _p[i][j] -= o._p[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  // exception throwing example
  for (auto i = 0; i < _rows; i++) {
    for (auto j = 0; j < _cols; j++) {
      _p[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& o) {
  if (_cols != o._rows)
    throw std::out_of_range(
        "the number of columns of the first matrix is not equal to the number "
        "of"
        " rows of the second matrix ");
  S21Matrix tmp = S21Matrix(this->_rows, o._cols);
  for (int i = 0; i < tmp._rows; i++) {
    for (int j = 0; j < tmp._cols; j++) {
      for (int v = 0; v < this->_cols; v++) {
        tmp._p[i][j] += this->_p[i][v] * o._p[v][j];
      }
    }
  }
  *this = tmp;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix tmp = S21Matrix(_cols, _rows);
  for (int i = 0; i < tmp._rows; i++) {
    for (int j = 0; j < tmp._cols; j++) {
      tmp._p[i][j] = _p[j][i];
    }
  }
  return tmp;
}

S21Matrix S21Matrix::looking_for_minor(int rows, int columns) {
  S21Matrix result = S21Matrix(_rows - 1, _cols - 1);
  int v = 0;
  for (int i = 0; i < _rows; i++) {
    if (i != rows) {
      int s = 0;
      for (int j = 0; j < _cols; j++) {
        if (j != columns) {
          result._p[v][s++] = _p[i][j];
        }
      }
      v++;
    }
  }
  return result;
}

double S21Matrix::Determinant() {
  if (_cols != _rows) {
    throw std::out_of_range("The matrix isn't square!");
  }
  double res = 0;
  if (_rows == 1) {
    res = _p[0][0];
  } else if (_rows == 2) {
    res = _p[0][0] * _p[1][1] - _p[0][1] * _p[1][0];
  } else {
    for (auto i = 0; i < _cols; i++) {
      S21Matrix minor_m = this->looking_for_minor(0, i);
      double tmp = minor_m.Determinant();
      double minor_2 = pow(-1, i) * _p[0][i] * tmp;
      res += minor_2;
    }
  }
  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  if (_cols != _rows) {
    throw std::out_of_range("The matrix isn't square!");
  }
  S21Matrix result = S21Matrix(_rows, _cols);
  for (int i = 0; i < result._rows; i++) {
    for (int j = 0; j < result._cols; j++) {
      S21Matrix minor = this->looking_for_minor(i, j);
      double deter = minor.Determinant();
      deter *= pow(-1, i + j);
      result._p[i][j] = deter;
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  S21Matrix result;
  double deter = this->Determinant();
  if (fabs(deter) <= 1e-7) {
    throw std::out_of_range("Matrix determinant is 0");
  } else {
    S21Matrix xxx_minor = this->CalcComplements();
    S21Matrix trans = xxx_minor.Transpose();
    trans.MulNumber(1 / deter);
    result.copy_obj(trans);
  }
  return result;
}

S21Matrix S21Matrix::operator+(const S21Matrix& o) {
  S21Matrix res(*this);
  res.SumMatrix(o);
  return res;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& o) {
  this->SumMatrix(o);
  return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix& o) {
  S21Matrix res(*this);
  res.SubMatrix(o);
  return res;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& o) {
  this->SubMatrix(o);
  return *this;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& o) {
  if (this->EqMatrix(o)) {
    return *this;
  }
  this->remove_matrix();
  this->_rows = o._rows;
  this->_cols = o._cols;
  this->create_matrix();
  this->copy_obj(o);
  return *this;
}

S21Matrix S21Matrix::operator*(const S21Matrix& o) {
  S21Matrix res(*this);
  res.MulMatrix(o);
  return res;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& o) {
  this->MulMatrix(o);
  return *this;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix res(*this);
  res.MulNumber(num);
  return res;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  this->MulNumber(num);
  return *this;
}

bool S21Matrix::operator==(const S21Matrix& o) {
  bool res = this->EqMatrix(o);
  return res;
}
double& S21Matrix::operator()(int row, int col) {
  if (row < 0 || col < 0 || row > _rows || col > _cols) {
    throw std::out_of_range("Index is outside the matrix");
  }
  return _p[row][col];
}

//  getters
int S21Matrix::get_rows() { return _rows; }
int S21Matrix::get_cols() { return _cols; }
double** S21Matrix::get_p() { return _p; }

//  setters
void S21Matrix::set_rows(int rows) {
  int index = 0;
  if (rows > this->_rows) {
    index = 1;
  }
  S21Matrix m1 = S21Matrix(*this);
  this->remove_matrix();
  this->_rows = rows;
  this->create_matrix();
  this->copy_matrix(m1, index);
}
void S21Matrix::set_cols(int cols) {
  int index = 0;
  if (cols > this->_cols) {
    index = 1;
  }
  S21Matrix m1 = S21Matrix(*this);
  this->remove_matrix();
  this->_cols = cols;
  this->create_matrix();
  this->copy_matrix(m1, index);
}
void S21Matrix::set_value(double value, int rows, int cols) {
  if (this->_rows >= rows || this->_cols >= cols) {
    _p[rows][cols] = value;
  }
}
void S21Matrix::set_from_0_to_future() {
  int k = 0;
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _p[i][j] = k++;
    }
  }
}
void S21Matrix::copy_obj(const S21Matrix& o) {
  this->_rows = o._rows;
  this->_cols = o._cols;
  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < this->_cols; j++) {
      this->_p[i][j] = o._p[i][j];
    }
  }
}

void S21Matrix::copy_matrix(const S21Matrix& o, int i) {
  if (i) {
    for (int i = 0; i < o._rows; i++) {
      for (int j = 0; j < o._cols; j++) {
        this->_p[i][j] = o._p[i][j];
      }
    }
  } else {
    for (int i = 0; i < this->_rows; i++) {
      for (int j = 0; j < this->_cols; j++) {
        this->_p[i][j] = o._p[i][j];
      }
    }
  }
}

void S21Matrix::create_matrix() {
  this->_p = new double*[this->_rows];
  for (int i = 0; i < this->_rows; i++) {
    this->_p[i] = new double[this->_cols]();
  }
}
void S21Matrix::remove_matrix() {
  if (this->_p) {
    for (auto i = 0; i < this->_rows; i++) {
      delete this->_p[i];
    }
    delete[] this->_p;
  }
}
void S21Matrix::m_out() {
  for (int i = 0; i < this->_rows; i++) {
    printf("\n");
    for (int j = 0; j < this->_cols; j++) {
      printf("%.2lf ", this->_p[i][j]);
    }
  }
  std::cout << "\n" << std::endl;
}
void S21Matrix::zero() {
  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < this->_cols; j++) {
      this->_p[i][j] = 0;
    }
  }
}
//
// int main() {
//    S21Matrix m1(3,3);
//    S21Matrix m2(3, 3);
//    m1.set_from_0_to_future();
//    m2.set_from_0_to_future();
//    m1.m_out();
//    m2.set_rows(4);
//    m2.m_out();
//    std::cout << m1.EqMatrix(m2);
//}