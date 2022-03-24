#ifndef PRIMITIVEVARS_HPP
#define PRIMITIVEVARS_HPP

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "geometry/vector.hpp"

#include "compressible.hpp"

class Compressible;

class PrimitiveVars {
public:
  double p;
  Vector2d u;
  double T;

  PrimitiveVars() {};
  PrimitiveVars(const double& P, const Vector2d& U, const double& t): p(P), u(U), T(t) {};
  ~PrimitiveVars() {};

  void zero() {p=0.; u = Vector2d(0., 0.); T=0.;}
  void one() {p=1.; u = Vector2d(1., 1.); T=1.;}

  static PrimitiveVars set(const Compressible& w);

  static const int nVars = 4;

  double& operator[](int k) {
    switch (k) {
    case 0: return p;
      break;
    case 1: return u.x;
      break;
    case 2: return u.y;
      break;
    case 3: return T;
      break;
    default: std::cout << "Out of range of primitive variables!" << std::endl;
      exit(10);
    }
  }

  double operator[](int k) const {
    switch (k) {
    case 0: return p;
      break;
    case 1: return u.x;
      break;
    case 2: return u.y;
      break;
    case 3: return T;
      break;
    default: std::cout << "Out of range of primitive variables!" << std::endl;
      exit(10);
    }
  }
};

inline PrimitiveVars operator+(const PrimitiveVars& a, const PrimitiveVars& b) {
  return PrimitiveVars(a.p+b.p, a.u+b.u, a.T+b.T);
}

inline PrimitiveVars operator-(const PrimitiveVars& a, const PrimitiveVars& b) {
  return PrimitiveVars(a.p-b.p, a.u-b.u, a.T-b.T);
}

template <typename S>
inline PrimitiveVars operator*(const PrimitiveVars& a, const S& b) {
  return PrimitiveVars(a.p*b, a.u*b, a.T*b);
}

template <typename S>
inline PrimitiveVars operator*(const S& b, const PrimitiveVars& a) {
  return PrimitiveVars(a.p*b, a.u*b, a.T*b);
}

inline PrimitiveVars operator*(const PrimitiveVars& a, const PrimitiveVars& b) {
  return PrimitiveVars(a.p*b.p, a.u*b.u, a.T*b.T);
}

template <typename S>
inline PrimitiveVars operator/(const PrimitiveVars& a, const S& b) {
  return PrimitiveVars(a.p/b, a.u/b, a.T/b);
}

inline PrimitiveVars operator+=(PrimitiveVars& a, const PrimitiveVars& b) {
  a.p+=b.p; a.u+=b.u; a.T+=b.T;
  return a;
}

inline PrimitiveVars operator-=(PrimitiveVars& a, const PrimitiveVars& b) {
  a.p-=b.p; a.u-=b.u; a.T-=b.T;
  return a;
}

template <typename S>
inline PrimitiveVars operator/=(PrimitiveVars& a, const S& b) {
  a.p/=b; a.u/=b; a.T/=b;
  return a;
}

#endif
