#include "setBoundaryConditions.hpp"

void setBoundaryConditions(CellField<Compressible>& w, const Grid& g, const Settings& setting,
			   const map<string, bcWithJacobian>& BC) {

  // nastaveni vstupni a vystupni hranice
#pragma omp parallel for
  for (int j=0; j<w.N(); j++) {
    // vstupni hranice
    const Face& fLeft = g.faceJ(0, j);

    auto it = BC.find(fLeft.name);   // auto -> map<string, bCondition>::const_iterator
    w[-1][j] = it->second.first(w[0][j], fLeft.s, setting);

    // vystupni hranice
    const Face& fRight = g.faceJ(w.M(), j);

    it = BC.find(fRight.name);
    w[w.M()][j] = it->second.first(w[w.M()-1][j], fRight.s, setting);
  }

  // nastveni dolni a horni hranice
#pragma omp parallel for
  for (int i=-w.gh(); i<w.M()+w.gh(); i++) {
    // dolni hranice
    const Face& fLower = g.faceI(i, 0);

    auto it = BC.find(fLower.name);
    if (it!=BC.end()) {
      w[i][-1] = it->second.first(w[i][0], fLower.s, setting);
    }
    else {
      cout << "Okrajova podminka nenalezena! Problem ve funkci setBoundarConditions." << endl;
      exit(14);
    }

    // horni hranice
    const Face& fUpper = g.faceI(i, w.N());

    it = BC.find(fUpper.name);
    w[i][w.N()] = it->second.first(w[i][w.N()-1], fUpper.s, setting);
  }
}
