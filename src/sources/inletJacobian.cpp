#include "inletJacobian.hpp"

Matrixd inletJacobian(const Compressible& wInside, const Vector2d& s, const Settings& setting) {
  double p = wInside.p();
  const double& p0 = setting.p0;
  const double& rho0 = setting.rho0;
  const double& alpha = setting.alpha;
  const double& kappa = Compressible::kappa;

  double ulx = wInside.rhoU.x / wInside.rho;
  double uly = wInside.rhoU.y / wInside.rho;

  double V = pow(p/p0, (1.-kappa)/kappa);

  double epsilon = 1.e-6;
  if(V<=1. + epsilon) V=1. + epsilon;

  double Ma = sqrt(2. / (kappa - 1.) * (V - 1.));
  double a = sqrt(kappa * p0 / rho0 * 1./V);
  
  double Dp[4];
  Dp[0] = (kappa-1.) * (pow(ulx, 2) + pow(uly, 2)) / 2.;
  Dp[1] = -(kappa-1.) * ulx;
  Dp[2] = -(kappa-1.) * uly;
  Dp[3] = kappa-1.;

  Matrixd J(Compressible::nVars);
  J.zero();

  double coeff = 1. / (a*a);
  // 1. radek matice
  J[0][0] = coeff * Dp[0];
  J[0][1] = coeff * Dp[1];
  J[0][2] = coeff * Dp[2];
  J[0][3] = coeff * Dp[3];

  coeff = 1./a * (Ma - V/Ma + Ma*(kappa-1.)/2.);
  // 2. radek matice
  J[1][0] = coeff *Dp[0] * cos(alpha);
  J[1][1] = coeff *Dp[1] * cos(alpha);
  J[1][2] = coeff *Dp[2] * cos(alpha);
  J[1][3] = coeff *Dp[3] * cos(alpha);
  
  // 3. radek matice
  J[2][0] = coeff *Dp[0] * sin(alpha);
  J[2][1] = coeff *Dp[1] * sin(alpha);
  J[2][2] = coeff *Dp[2] * sin(alpha);
  J[2][3] = coeff *Dp[3] * sin(alpha);

  coeff = (V - kappa + 1.) / (kappa - 1.);
  // 4. radek matice
  J[3][0] = coeff * Dp[0];
  J[3][1] = coeff * Dp[1];
  J[3][2] = coeff * Dp[2];
  J[3][3] = coeff * Dp[3];
 
  return J;
}
