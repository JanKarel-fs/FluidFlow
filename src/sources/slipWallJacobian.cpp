#include "slipWallJacobian.hpp"

Matrixd slipWallJacobian(const Compressible& wInside, const Vector2d& s, const Settings& setting) {
  Matrixd J(Compressible::nVars);

  J.ones();

  Vector2d n = s / s.length();

  J[1][1] = 1. - 2. * n.x*n.x;   J[1][2] = -2. * n.x*n.y;
  J[2][1] = -2. * n.x*n.y;       J[2][2] = 1. - 2. * n.y*n.y;

  return J;
}
