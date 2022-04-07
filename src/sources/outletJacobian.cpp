#include "outletJacobian.hpp"

Matrixd outletJacobian(const Compressible& wInside, const Vector2d& s, const Settings& setting) {
  const double& rho = wInside.rho;
  const Vector2d& rhoU = wInside.rhoU;
  
  Matrixd J(Compressible::nVars);

  J.zero();

  J[0][0] = 1.;
  J[1][1] = 1.;
  J[2][2] = 1.;

  J[3][0] = -0.5 * (pow(rhoU.x, 2) + pow(rhoU.y, 2)) / (pow(rho, 2));
  J[3][1] = rhoU.x / rho;
  J[3][2] = rhoU.y / rho;

  return J;
}
