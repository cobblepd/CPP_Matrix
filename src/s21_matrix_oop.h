#ifndef s21_MATRIX_OOP_H_
#define S21_MATRIX_OOP_H_
#define DEFAULT_ROWS 1
#define DEFAULT_COLUMNS 1
using namespace std;

class S21Matrix {
 private:
  int rows_ = 0, cols_ = 0;
  double **matrix_ = nullptr;

  void NullRow(int row_i, int cols, int start = 0) noexcept;
  void Alloc(int rows, int cols, bool set_null);
  void FreeMatrix() noexcept;
  void SumSubMatrix(const S21Matrix &other, bool is_sum) noexcept;

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other);
  ~S21Matrix() noexcept;

  int GetRows() const noexcept;
  int GetCols() const noexcept;
  void SetRows(int rows);
  void SetCols(int cols);

  // operations
  bool EqMatrix(const S21Matrix &other) noexcept;
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num) noexcept;
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  double Minor(int i, int j);
  S21Matrix InverseMatrix();

  // operator overloads
  double &operator()(int i, int j);
  double operator()(int i, int j) const;
  S21Matrix &operator=(const S21Matrix &other);
  S21Matrix &operator=(S21Matrix &&other);
  S21Matrix operator+(const S21Matrix &other) const;
  S21Matrix operator+() const noexcept;
  S21Matrix operator-(const S21Matrix &other) const;
  S21Matrix operator-() const;
  S21Matrix operator*(const S21Matrix &other) const;
  S21Matrix operator*(const double num) const;
  friend S21Matrix operator*(double num, const S21Matrix &other);
  bool operator==(const S21Matrix &other) noexcept;
  S21Matrix operator+=(const S21Matrix &other);
  S21Matrix operator-=(const S21Matrix &other);
  S21Matrix operator*=(const S21Matrix &other);
  S21Matrix operator*=(const double &num) noexcept;
};

#endif
