#include "storeNormResidue.hpp"

void storeNormResidue(const CellField<Compressible>& res, const Grid& g, const int& it) {
  Compressible L2(0., Vector2d(0., 0.), 0);
  Compressible Linf(0., Vector2d(0., 0.), 0);

  for (int i=0; i<res.M(); i++) {
    for (int j=0; j<res.N(); j++) {
      const Compressible& R = res[i][j];
      const double& volume = g.volume(i, j);

      Linf = Compressible::max(Compressible::fabs(R), Linf);

      L2 += R * R * volume;
    }
  }

  L2 = Compressible::sqrt(L2);

  L2.rho = log10(L2.rho);
  L2.rhoU.x = log10(L2.rhoU.x);
  L2.rhoU.y = log10(L2.rhoU.y);
  L2.e = log10(L2.e);

  Linf.rho = log10(Linf.rho);
  Linf.rhoU.x = log10(Linf.rhoU.x);
  Linf.rhoU.y = log10(Linf.rhoU.y);
  Linf.e = log10(Linf.e);

  ofstream fout("results/L2.txt", ios::app);
  fout << it << " " << L2.rho << " " << L2.rhoU.x << " " << L2.rhoU.y << " " << L2.e << endl;
  fout.close();

  fout.open("results/Linf.txt", ios::app);
  fout << it << " " << Linf.rho << " " << Linf.rhoU.x << " " << Linf.rhoU.y << " " <<Linf.e << endl;
  fout.close();
}
