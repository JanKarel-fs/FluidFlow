#include "saveResults.hpp"

void saveResults(const CellField<Compressible>& w, const Grid& g) {
  NodeField<Compressible> wnd(g);

  for (int i=0; i<wnd.M(); i++) {
    for (int j=0; j<wnd.N(); j++) {
      double volume = g.volume(i, j) + g.volume(i-1, j) + g.volume(i-1, j-1) + g.volume(i, j-1);

      wnd[i][j] = w[i][j] * g.volume(i, j) + w[i-1][j] * g.volume(i-1, j);
      wnd[i][j] += w[i-1][j-1] * g.volume(i-1, j-1) + w[i][j-1] * g.volume(i, j-1);

      wnd[i][j] /= volume;
    }
  }

  ofstream fout("results/results.dat");
  fout << "VARIABLES = \"x\", \"y\", \"rho\", \"u\", \"v\", \"e\", \"p\", \"a\", \"Ma\"" << endl;
  fout << "ZONE T = \"node_values\", I = " << wnd.M() << ", J = " << wnd.N() << ", F = POINT" << endl;

  for (int j=0; j<wnd.N(); j++) {
    for (int i=0; i<wnd.M(); i++) {
      fout << g.x(i, j) << " " << g.y(i, j) << " " << wnd[i][j].rho;
      fout << " " << wnd[i][j].rhoU.x/wnd[i][j].rho << " " << wnd[i][j].rhoU.y/wnd[i][j].rho;
      fout << " " << wnd[i][j].e << " " << wnd[i][j].p() << " " << wnd[i][j].a();
      fout << " " << wnd[i][j].Ma() << endl;
    }
  }

  fout.close();
}
