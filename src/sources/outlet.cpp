#include "outlet.hpp"

Compressible outlet(const Compressible& wInside, const Vector2d& s, const Settings& setting) {
  const double& p0 = setting.p0;
  const double& Ma2is = setting.Ma2is;
  const double& kappa = Compressible::kappa;

  const double& rho = wInside.rho;
  const Vector2d& rhoU = wInside.rhoU;

  double p = p0 * pow(1. + (kappa - 1.) / 2. * pow(Ma2is, 2), kappa / (1. - kappa));

  double e = p / (kappa - 1.) + 0.5 * (pow(rhoU.x, 2) + pow(rhoU.y, 2)) / rho;

  return Compressible(rho, rhoU, e);
}
