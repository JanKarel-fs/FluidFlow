#include "primitiveVars.hpp"

PrimitiveVars PrimitiveVars::set(const Compressible& w) {
  double P = w.p();
  Vector2d U = w.rhoU / w.rho;
  double t = w.T();

  return PrimitiveVars(P, U, t);
}
