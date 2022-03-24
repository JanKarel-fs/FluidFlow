#ifndef COMPUTERESIDUE_HPP
#define COMPUTERESIDUE_HPP

#include <iostream>
#include <cstdlib>
#include "grid.hpp"
#include "cellField.hpp"
#include "grad.hpp"
#include "limiter.hpp"
#include "computeResidueConv.hpp"
#include "computeResidueDiss.hpp"
#include "../geometry/vector.hpp"
#include "sources/settings.hpp"
#include <omp.h>

using namespace std;

template <typename var>
void computeResidue(const CellField<var>& w, const Grid& g, CellField<var>& res,
		    const Settings& setting) {
  // vynulovani residua
#pragma omp parallel for
  for (int i=0; i<res.M(); i++) {
    for (int j=0; j<res.N(); j++) {
      res[i][j].zero();
    }
  }

  switch (setting.convection) {
  case 0: break;      // konvekce je vypnuta nedelame nic
  case 1: computeResidueConv(w, g, res);
    break;
  default: cout << "No such possibility for convection!" << endl;
    exit(15);
  }

  switch (setting.diffusion) {
  case 0: break;      // difuze je vypnuta nedelame nic
  case 1: computeResidueDiss(w, g, res);
    break;
  default: cout << "No such possibility for convection!" << endl;
    exit(15);
  }

  // deleni objemem bunky
  for (int i=0; i<res.M(); i++) {
    for (int j=0; j<res.N(); j++) {
      res[i][j] /= g.volume(i, j);
    }
  }

}

#endif
