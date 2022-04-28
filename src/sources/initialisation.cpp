#include "initialisation.hpp"

void initialisation(CellField<Compressible>& w, const Grid& g, const Settings& setting) {

  system("mkdir results");  // Linux i Windows
  // system("md results")  // Windows
  system("rm results/*");  // Linux
  // system("del results\*"); // Windows (s otaznikem, lomitko vs. zpetne lomitko)

  switch (setting.solver) {
  case 1: {
    step<Compressible> = stepExplicit<Compressible>;
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
  }
    break;

  case 2: {
    step<Compressible> = stepImplicit<Compressible>;
    switch (setting.flux) {
    case 1:
      Compressible::fluxImplicit = Compressible::UpwindImplicit;
      break;
    case 2:
      Compressible::fluxImplicit = Compressible::RusanovImplicit;
      break;
    default:
      cout << "Takovy tok neexistuje!" << endl;
      exit(0);
    }

    switch (setting.temporalOrder) {
    case 1: timeIncrement<Compressible> = timeIncrementFirstOrder<Compressible>;
      break;
    case 2:  timeIncrement<Compressible> = timeIncrementSecondOrder<Compressible>;
      break;
    default:
      cout << "No such temporal order!" << endl;
      exit(0);
    }
  }
    break;

  default:
    cout << "There is no possible choice for a solver!" << endl;
    cout << "Possibilities are 1 for explicit solver, 2 for implicit solver." << endl;
    exit(1);
  }

  switch (setting.spatialOrder) {
  case 1:
    limiter<Compressible> = zeroLimiter<Compressible>;
    grad<Compressible> = zeroGrad<Compressible>;
    break;
  case 2:
    switch (setting.limiter) {
    case 1:
      limiter<Compressible> = barthJespersen<Compressible>;
      break;
    case 2:
      limiter<Compressible> = venkatakrishnan<Compressible>;
      break;
    default:
      cout << "No such possibility for limiter" << endl;
      exit(1);
    }
    grad<Compressible> = gradLSM<Compressible>;
    break;
  default:
    cout << "No such spatial order!" << endl;
    exit(1);
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
