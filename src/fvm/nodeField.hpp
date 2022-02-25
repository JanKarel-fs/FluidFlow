#ifndef NODEFIELD_HPP
#define NODEFIELD_HPP

#include <iostream>
#include <cstdlib>
#include "../geometry/field.hpp"
#include "grid.hpp"

using namespace std;

template <typename T>
class NodeField {
  Field2<T> data;

  int m;
  int n;
  int ghost;

public:
  NodeField() {};

  NodeField(const int M, const int N, const int gh): m(M), n(N), ghost(gh) {
    data.allocate(-ghost, m+ghost, -ghost, n+ghost);
  }

  NodeField(const Grid& g) {
    m = g.Mnd();
    n = g.Nnd();
    ghost = g.ghostNum();

    data.allocate(-ghost, m+ghost, -ghost, n+ghost);
  };

  ~NodeField() {};

  inline T* operator[](const int i) const {
    return data[i];
  }

  int M() const {return m;}
  int N() const {return n;}
  int gh() const {return ghost;}
  int Imin() const {return data.Imin();}
  int Imax() const {return data.Imax();}
  int Jmin() const {return data.Jmin();}
  int Jmax() const {return data.Jmax();}

  NodeField& operator=(const NodeField& cf) {
    if (ghost!=cf.gh() || m!=cf.M() || n!=cf.N()) {
      cout << "Nelze priradit hodnoty z pole NodeField jine velikosti!" << endl;
      exit(2);
    }

    for (int i=cf.Imin(); i<cf.Imax(); i++) {
      for (int j=cf.Jmin(); j<cf.Jmax(); j++) {
	data[i][j] = cf[i][j];
      }
    }

    return *this;
  }
};

template <typename T>
inline NodeField<T> operator+(const NodeField<T>& a, const NodeField<T>& b) {
  if (a.gh()!=b.gh() || a.M()!=b.M() || a.N()!=b.N()) {
    cout << "Nelze scitat pole NodeField ruznych velikosti!" << endl;
    exit(3);
  }

  NodeField<T> c(a.M(), a.N(), a.gh());

  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      c[i][j] = a[i][j] + b[i][j];
    }
  }

  return c;
}

template <typename T>
inline NodeField<T> operator-(const NodeField<T>& a, const NodeField<T>& b) {
  if (a.gh()!=b.gh() || a.M()!=b.M() || a.N()!=b.N()) {
    cout << "Nelze odecitat pole NodeField ruznych velikosti!" << endl;
    exit(3);
  }

  NodeField<T> c(a.M(), a.N(), a.gh());

  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      c[i][j] = a[i][j] - b[i][j];
    }
  }

  return c;
}

template <typename T, typename S>
inline NodeField<T> operator*(const NodeField<T>& a, const S& b) {
  NodeField<T> c(a.M(), a.N(), a.gh());

  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      c[i][j] = a[i][j] * b;
    }
  }

  return c;
}


template <typename T, typename S>
inline NodeField<T> operator*(const S& b, const NodeField<T>& a) {
  NodeField<T> c(a.M(), a.N(), a.gh());

  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      c[i][j] = a[i][j] * b;
    }
  }

  return c;
}

template <typename T, typename S>
inline NodeField<T> operator/(const NodeField<T>& a, const S& b) {
  NodeField<T> c(a.M(), a.N(), a.gh());

  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      c[i][j] = a[i][j] / b;
    }
  }

  return c;
}

#endif
