#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cassert>

#include "vector.hpp"
#include <iostream>

using namespace std;

template <typename T>
class Matrix {
public:
  int M;
  int N;

  T *field;
    
  Matrix(): M(0), N(0) {allocated = false;};
  Matrix(int i): M(i), N(i) {field = new T[M*N]; allocated = true;};
  Matrix(int i, int j): M(i), N(j) {field = new T[M*N]; allocated = true;};
  ~Matrix() {if (allocated) delete [] field;};

  Matrix(const Matrix& A) {
    resize(A.M, A.N);

    int totalSize = M*N;

    for (int i=0; i<totalSize; i++) {
      field[i] = A.field[i];
    }
  };

  void resize(const int i, const int j) {
    if (allocated) {
      delete [] field;
    }

    M = i;
    N = j;

    field = new T[M*N];
    allocated = true;
  }

  Matrix& operator=(const Matrix& A) {
    resize(A.M, A.N);

    int totalSize = M*N;

    for (int i=0; i<totalSize; i++) {
      field[i] = A.field[i];
    }
    return *this;
  }
    
  void zero();
  void ones();
  void print();
  Matrix transpose();

  inline T* operator[](const int i) {assert (0 <= i && i < M); return field + i*N;}
  inline const T* operator[](const int i) const {assert (0 <= i && i < M); return field + i*N;}

private:
  bool allocated;
};

typedef Matrix<double> Matrixd;

template <typename T>
inline Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b) {    
  Matrix<T> c(a.M, a.N);
  int size = a.M * a.N;
    
  if ((a.M == b.M) && (a.N == b.N)) {
    for (int i=0; i<size; i++) {
      c.field[i] = a.field[i] + b.field[i];
    }
  }
  else {
    cout << "Cannot " << endl;
    exit(0);
  }
    
  return c;
}

template <typename T>
inline Matrix<T> operator-(const Matrix<T>& a, const Matrix<T>& b) {    
  Matrix<T> c(a.M, a.N);
  int size = a.M * a.N;
    
  if ((a.M == b.M) && (a.N == b.N)) {
    for (int i=0; i<size; i++) {
      c.field[i] = a.field[i] - b.field[i];
    }
  }
  else {
    cout << "Cannot " << endl;
    exit(0);
  }
    
  return c;
}

template <typename T>
inline Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b) {    
  Matrix<T> c(a.M, b.N);
  c.zero();
    
  if (a.N == b.M) {
    for (int i=0; i<a.M; i++) {
      for (int j=0; j<b.N; j++) {
	for (int k=0; k<a.N; k++) {
	  c.field[i*b.N + j] += a.field[i*a.N + k] * b.field[k*b.N + j];
	}
      }
    }
  }
  else {
    cout << "Cannot " << endl;
    exit(0);
  }
    
  return c;
}

template <typename T, typename S>
inline Matrix<T> operator*(const Matrix<T>& a, S b) {    
  Matrix<T> c(a.M, a.N);
  int size = a.M * a.N;
    
  for (int i=0; i<size; i++) {
    c.field[i] = a.field[i] * b;
  }
    
  return c;
}

template <typename T, typename S>
inline Matrix<T> operator*(S b, const Matrix<T>& a) {    
  Matrix<T> c(a.M, a.N);
  int size = a.M * a.N;
    
  for (int i=0; i<size; i++) {
    c.field[i] = a.field[i] * b;
  }
    
  return c;
}

template <typename T, typename S>
inline Matrix<T> operator/(const Matrix<T>& a, S b) {
  Matrix<T> c(a.M, a.N);
  int size = a.M * a.N;
    
  for (int i=0; i<size; i++) {
    c.field[i] = a.field[i] / b;
  }
    
  return c;
}

template <typename T>
void Matrix<T>::zero() {
  int size = M * N;
  
  for (int i=0; i<size; i++) {
    field[i] = 0;
  }
}

template <typename T>
void Matrix<T>::ones() {
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      if (i == j) {
	field[i*N + j] = 1;
      }
      else {
	field[i*N + j] = 0;
      }
    }
  }
}

template <typename T>
Matrix<T> Matrix<T>::transpose() {    
  Matrix<T> c(N, M);
    
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      c.field[j*M + i] = field[i*N + j];
    }
  }
    
  return c;
}

template <typename T>
inline Vector2<T> operator*(const Matrix<T>& A, const Vector2<T>& b) {

  if (A.M != 2 || A.N != 2) {
    cout << "Cannot multilply a matrix " << A.M << "x" << A.N << " by a Vector3!" << endl;
    cout << "Matrix has to be 2x2!" << endl;
    exit(0);
  }
  
  T c[2], d[2];
  c[0] = b.x;  c[1] = b.y;
  
  for (int i=0; i<2; i++) {
    d[i] = 0;
    for (int j=0; j<2; j++) {
      d[i] += A[i][j] * c[j];
    }
  }

  return Vector2<T>(d[0], d[1]);
}

template <typename T>
void Matrix<T>::print() {
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      cout << field[i*N + j] << "\t";
    }
    cout << endl;
  }
}

#endif
