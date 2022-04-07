#include "noSlipWallJacobian.hpp"

Matrixd noSlipWallJacobian(const Compressible& wInside, const Vector2d& s, const Settings& setting) {
  Matrixd J(Compressible::nVars);

  J.ones();

  J[1][1] = -1.;
  J[2][2] = -1.;

  return J;
}
