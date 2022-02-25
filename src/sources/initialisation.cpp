#include "initialisation.hpp"

void initialisation(CellField<Compressible>& w, const Grid& g, const Settings& setting) {

  system("mkdir results");  // Linux i Windows
  // system("md results")  // Windows
  system("rm results/*");  // Linux
  // system("del results\*"); // Windows (s otaznikem, lomitko vs. zpetne lomitko)

 switch (setting.flux) {
  case 1:
    Compressible::flux = Compressible::Upwind;
    break;
  case 2:
    Compressible::flux = Compressible::Rusanov;
    break;
  default:
    cout << "Takovy tok neexistuje!" << endl;
    exit(0);
  }

 Compressible::kappa = setting.kappa;
 Compressible::R = setting.R;
 Compressible::Pr = setting.Pr;
 Compressible::cv = Compressible::R / (Compressible::kappa - 1.);
 Compressible::cp = Compressible::cv * Compressible::kappa;

 const double& rho = setting.rhoInit;
 const Vector2d& u = setting.uInit;
 double e = setting.pInit / (setting.kappa - 1.) + 0.5 * rho *(pow(u.x, 2) + pow(u.y, 2));

  for (int i=0; i<w.M(); i++) {
    for (int j=0; j<w.N(); j++) {
      w[i][j] = Compressible(rho, rho*u, e);
    }
  }
}