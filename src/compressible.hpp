#ifndef COMPRESSIBLE_HPP
#define COMPRESSIBLE_HPP

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "geometry/vector.hpp"

class Compressible {
public:
  double rho;
  Vector2d rhoU;
  double e;

  Compressible() {};
  Compressible(const double& Rho, const Vector2d& RhoU, const double& E): rho(Rho), rhoU(RhoU), e(E) {};
  ~Compressible() {};

  double p() const;
  double a() const;
  double Ma() const;

  void zero() {rho=0.; rhoU = Vector2d(0., 0.); e=0.;}
  void one() {rho=1.; rhoU = Vector2d(1., 1.); e=1.;}

  static Compressible (*flux)(const Compressible& wl, const Compressible& wr, const Vector2d& s);
  static Compressible Upwind(const Compressible& wl, const Compressible& wr, const Vector2d& s);
  static Compressible Rusanov(const Compressible& wl, const Compressible& wr, const Vector2d& s);
  static Compressible fabs(const Compressible& w);
  static Compressible sqrt(const Compressible& w);
  static Compressible max(const Compressible& a, const Compressible& b);
  static Compressible min(const Compressible& a, const Compressible& b);
  
  static double kappa;
  static double cp;
  static double cv;
  static double R;
  static double Pr;
  static const int nVars = 4;

  double& operator[](int k) {
    switch (k) {
    case 0: return rho;
      break;
    case 1: return rhoU.x;
      break;
    case 2: return rhoU.y;
      break;
    case 3: return e;
      break;
    default: std::cout << "Out of range of compressible variables!" << std::endl;
      exit(10);
    }
  }

  double operator[](int k) const {
    switch (k) {
    case 0: return rho;
      break;
    case 1: return rhoU.x;
      break;
    case 2: return rhoU.y;
      break;
    case 3: return e;
      break;
    default: std::cout << "Out of range of compressible variables!" << std::endl;
      exit(10);
    }
  }
};

inline Compressible operator+(const Compressible& a, const Compressible& b) {
  return Compressible(a.rho+b.rho, a.rhoU+b.rhoU, a.e+b.e);
}

inline Compressible operator-(const Compressible& a, const Compressible& b) {
  return Compressible(a.rho-b.rho, a.rhoU-b.rhoU, a.e-b.e);
}

template <typename S>
inline Compressible operator*(const Compressible& a, const S& b) {
  return Compressible(a.rho*b, a.rhoU*b, a.e*b);
}

template <typename S>
inline Compressible operator*(const S& b, const Compressible& a) {
  return Compressible(a.rho*b, a.rhoU*b, a.e*b);
}

inline Compressible operator*(const Compressible& a, const Compressible& b) {
  return Compressible(a.rho*b.rho, a.rhoU*b.rhoU, a.e*b.e);
}

template <typename S>
inline Compressible operator/(const Compressible& a, const S& b) {
  return Compressible(a.rho/b, a.rhoU/b, a.e/b);
}

inline Compressible operator+=(Compressible& a, const Compressible& b) {
  a.rho+=b.rho; a.rhoU+=b.rhoU; a.e+=b.e;
  return a;
}

inline Compressible operator-=(Compressible& a, const Compressible& b) {
  a.rho-=b.rho; a.rhoU-=b.rhoU; a.e-=b.e;
  return a;
}

template <typename S>
inline Compressible operator/=(Compressible& a, const S& b) {
  a.rho/=b; a.rhoU/=b; a.e/=b;
  return a;
}

#endif
