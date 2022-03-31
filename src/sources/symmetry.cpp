#include "symmetry.hpp"

Compressible symmetry(const Compressible& wInside, const Vector2d& s, const Settings& setting) {
  const double& rho = wInside.rho;
  const double& e = wInside.e;

  Vector2d n = s / s.length();

  Vector2d rhoU = wInside.rhoU - 2. * dot(wInside.rhoU, n) * n;

  return Compressible(rho, rhoU, e);
}
