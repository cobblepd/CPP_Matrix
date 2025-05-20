#include "s21_matrix_oop.h"

#include <new>
#include <stdexcept>

// ---------------------  CONSTRUCTORS  ---------------------

S21Matrix::S21Matrix() {
  rows_ = DEFAULT_ROWS;
  cols_ = DEFAULT_COLUMNS;
  Alloc(DEFAULT_ROWS, DEFAULT_COLUMNS, true);
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows > 0 && cols > 0) {
    rows_ = rows;
    cols_ = cols;
    Alloc(rows, cols, true);
  } else {
    throw invalid_argument("Row count or Column count is less than 1");
  }
}

S21Matrix::S21Matrix(const S21Matrix &other) { *this = other; }

S21Matrix::S21Matrix(S21Matrix &&other) { *this = std::move(other); }

void S21Matrix::Alloc(int rows, int cols, bool set_null) {
  matrix_ = new double *[rows];
  for (int i = 0; i < rows; ++i) {
    matrix_[i] = nullptr;
  }
  try {
    for (int i = 0; i < rows; ++i) {
      matrix_[i] = new double[cols];
      if (set_null) NullRow(i, cols);
    }
  } catch (const std::bad_alloc &e) {
    FreeMatrix();
    throw e;
  }
}

void S21Matrix::FreeMatrix() noexcept {
  for (int i = 0; i < rows_ && matrix_; ++i) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
}

S21Matrix::~S21Matrix() noexcept { FreeMatrix(); }

// ---------------------  GETTERS/SETTERS  ---------------------

void S21Matrix::SetCols(int cols) {
  if (cols > 0) {
    int oldCols = cols_;
    cols_ = cols;
    for (int i = 0; i < rows_; ++i) {
      double *newRow = new double[cols];
      for (int j = 0; j < oldCols && j < cols; ++j) {
        newRow[j] = matrix_[i][j];
      }
      delete[] matrix_[i];
      matrix_[i] = newRow;
      NullRow(i, cols, oldCols);
    }
  } else {
    throw invalid_argument("Column count is less than 1");
  }
}

void S21Matrix::SetRows(int rows) {
  if (rows > 0) {
    double **newMatrix = new double *[rows];
    for (int i = 0; i < rows_ && i < rows; ++i) {
      newMatrix[i] = matrix_[i];
    }
    if (rows > rows_) {
      delete[] matrix_;
      matrix_ = newMatrix;
      for (int i = rows_; i < rows; ++i) {
        matrix_[i] = new double[cols_];
        NullRow(i, cols_);
        rows_ = i + 1;
      }
    } else {
      for (int i = rows; i < rows_; ++i) {
        delete[] matrix_[i];
      }
      delete[] matrix_;
      matrix_ = newMatrix;
      rows_ = rows;
    }
  } else {
    throw invalid_argument("Row count is less than 1");
  }
}

int S21Matrix::GetRows() const noexcept { return rows_; };
int S21Matrix::GetCols() const noexcept { return cols_; };

// ---------------------  OPERATIONS  ---------------------

void S21Matrix::NullRow(int row_i, int cols, int start) noexcept {
  for (int j = start; j < cols; ++j) {
    matrix_[row_i][j] = 0;
  }
}

bool S21Matrix::EqMatrix(const S21Matrix &other) noexcept {
  bool result = true;

  if (rows_ != other.rows_ || cols_ != other.cols_) {
    result = false;
  }
  if (result) {
    for (int i = 0; i < rows_ && result; ++i) {
      for (int j = 0; j < cols_ && result; ++j) {
        result = abs(other.matrix_[i][j] - matrix_[i][j]) <= 1E-6;
      }
    }
  }

  return result;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    SumSubMatrix(other, true);
  } else {
    throw invalid_argument("Matrices have different sizes in Sum");
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    SumSubMatrix(other, false);
  } else {
    throw invalid_argument("Matrices have different sizes in Sub");
  }
}

void S21Matrix::SumSubMatrix(const S21Matrix &other, bool is_sum) noexcept {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (is_sum) {
        matrix_[i][j] += other.matrix_[i][j];
      } else {
        matrix_[i][j] -= other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::MulNumber(const double num) noexcept {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ == other.rows_) {
    S21Matrix temp = S21Matrix(*this);
    SetCols(other.cols_);
    for (int i = 0; i < temp.rows_; ++i) {
      for (int j = 0; j < other.cols_; ++j) {
        double s = 0;
        for (int k = 0; k < temp.cols_; ++k) {
          s += temp.matrix_[i][k] * other.matrix_[k][j];
        }
        matrix_[i][j] = s;
      }
    }
  } else {
    throw invalid_argument(
        "Matrix A's Column count differs from Matrix B's Row count");
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix transposed = S21Matrix(cols_, rows_);
  for (int i = 0; i < cols_; ++i) {
    for (int j = 0; j < rows_; ++j) {
      transposed.matrix_[i][j] = matrix_[j][i];
    }
  }

  return transposed;
}

S21Matrix S21Matrix::CalcComplements() {
  S21Matrix compMatrix(rows_, cols_);
  if (cols_ == rows_ && cols_ > 1) {
    int sign = 1;
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        compMatrix.matrix_[i][j] = sign * Minor(i, j);
        sign *= -1;
      }
    }
  } else if (cols_ != rows_) {
    throw logic_error("Matrix is not a Square Matrix");
  } else {
    throw logic_error("Can't find 1x1 matrix's Minor");
  }

  return compMatrix;
}

double S21Matrix::Determinant() {
  double determinant = 0;

  if (cols_ == rows_ && rows_ > 1) {
    S21Matrix triang = S21Matrix(*this);
    bool iszero = false;

    for (int p = 0; p < rows_ - 1 && !iszero; ++p) {  // p means pivot
      for (int i = p + 1; i < rows_ && !iszero; ++i) {
        double m = -triang.matrix_[i][p] / triang.matrix_[p][p];
        iszero = triang.matrix_[p][p] == 0;
        for (int j = 0; j < cols_ && !iszero; ++j) {
          triang.matrix_[i][j] += triang.matrix_[p][j] * m;
        }
      }
    }

    determinant = iszero ? 0 : triang.matrix_[0][0];
    for (int i = 1; i < rows_ && !iszero; ++i) {
      determinant *= triang.matrix_[i][i];
    }
  } else if (rows_ == 1) {
    determinant = matrix_[0][0];
  } else {
    throw logic_error("Matrix is not a Square Matrix");
  }

  return determinant;
}

double S21Matrix::Minor(int i, int j) {
  S21Matrix minorm(rows_ - 1, cols_ - 1);
  for (int oi = 0, mi = 0; oi < rows_; ++oi) {
    for (int oj = 0, mj = 0; oj < cols_; ++oj) {
      mi = oi < i ? oi : oi - 1;
      mj = oj < j ? oj : oj - 1;
      if (oi != i && oj != j) {
        minorm.matrix_[mi][mj] = matrix_[oi][oj];
      }
    }
  }

  return minorm.Determinant();
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = Determinant();
  S21Matrix inverse(rows_, cols_);

  if (det != 0 && rows_ > 1) {
    inverse = CalcComplements().Transpose();
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        inverse.matrix_[i][j] /= det;
      }
    }
  } else if (rows_ == 1) {
    inverse.matrix_[0][0] = 1.0 / matrix_[0][0];
  } else {
    throw logic_error("Determinant of the Matrix is zero");
  }

  return inverse;
}

// ---------------------  OPERATOR OVERLOADS  ---------------------

double &S21Matrix::operator()(int i, int j) {
  double *rtn = nullptr;
  if (i >= 0 && i < rows_ && j >= 0 && j < cols_) {
    rtn = &matrix_[i][j];
  } else {
    throw invalid_argument("Indexes are out of the scope of this Matrix");
  }

  return *rtn;
}

double S21Matrix::operator()(int i, int j) const {
  double rtn = 0;
  if (i >= 0 && i < rows_ && j >= 0 && j < cols_) {
    rtn = matrix_[i][j];
  } else {
    throw invalid_argument("Indexes are out of the scope of this Matrix");
  }

  return rtn;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (this != &other && other.rows_ > 0 && other.cols_ > 0) {
    SetRows(other.rows_);
    SetCols(other.cols_);
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  } else if (other.rows_ == 0 || other.cols_ == 0) {
    throw invalid_argument("Can't copy empty array");
  }
  return *this;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) {
  if (this != &other && other.rows_ > 0 && other.cols_ > 0) {
    FreeMatrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;
    other.matrix_ = nullptr;
    other.rows_ = 0;
    other.cols_ = 0;
  } else if (other.rows_ == 0 || other.cols_ == 0) {
    throw invalid_argument("Can't move empty array");
  }
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) const {
  S21Matrix sum(*this);
  sum.SumMatrix(other);
  return sum;
}

S21Matrix S21Matrix::operator+() const noexcept { return *this; };

S21Matrix S21Matrix::operator-(const S21Matrix &other) const {
  S21Matrix sub(*this);
  sub.SubMatrix(other);
  return sub;
}

S21Matrix S21Matrix::operator-() const { return -1 * (*this); };

S21Matrix S21Matrix::operator*(const S21Matrix &other) const {
  S21Matrix mul(*this);
  mul.MulMatrix(other);
  return mul;
}

S21Matrix S21Matrix::operator*(const double num) const {
  S21Matrix mul(*this);
  mul.MulNumber(num);
  return mul;
}

S21Matrix operator*(double num, const S21Matrix &other) { return other * num; };

bool S21Matrix::operator==(const S21Matrix &other) noexcept {
  return EqMatrix(other);
};

S21Matrix S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const double &num) noexcept {
  MulNumber(num);
  return *this;
}
