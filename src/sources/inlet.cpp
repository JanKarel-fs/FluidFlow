#include "inlet.hpp"

Compressible inlet(const Compressible& wInside, const Vector2d& s, const Settings& setting) {
  double p = wInside.p();
  const double& p0 = setting.p0;
  const double& rho0 = setting.rho0;
  const double& alpha = setting.alpha;
  const double& kappa = Compressible::kappa;

  double Ma = sqrt(2. / (kappa - 1.) * (pow(p/p0, (1.-kappa)/kappa) - 1.));
  double rho = rho0 * pow(p/p0, 1./kappa);
  double a = sqrt(kappa * p / rho);
  double rhoU_vel = rho * Ma * a;

  Vector2d rhoU(rhoU_vel * cos(alpha), rhoU_vel * sin(alpha));

  double e = p / (kappa - 1.) + 0.5 * (pow(rhoU.x, 2) + pow(rhoU.y, 2)) / rho;

  return Compressible(rho, rhoU, e);
}
