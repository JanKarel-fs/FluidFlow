#ifndef COMPUTERESIDUEIMPLICIT_HPP
#define COMPUTERESIDUEIMPLICIT_HPP

#include <iostream>
#include <cstdlib>
#include "grid.hpp"
#include "cellField.hpp"
#include "grad.hpp"
#include "limiter.hpp"
#include "computeResidueConvImplicit.hpp"
#include "computeResidueDissImplicit.hpp"
#include "../geometry/vector.hpp"
#include "sources/settings.hpp"
#include "sources/typedefs.hpp"
#include "sources/linearSolver.hpp"
#include <omp.h>

using namespace std;

template <typename var>
void computeResidueImplicit(const CellField<var>& w, const Grid& g,
			    const map<string, bcWithJaconian>& BC, const double& dt,
			    LinearSolver<var>& linSolver, const Settings& setting) {
  timeIncrement<var>(linSolver, dt);

  switch (setting.convection) {
  case 0: break;      // konvekce je vypnuta nedelame nic
  case 1: computeResidueConvImplicit(w, g, BC, linSolver, setting);
    break;
  default: cout << "No such possibility for convection!" << endl;
    exit(15);
  }

  switch (setting.diffusion) {
  case 0: break;      // difuze je vypnuta nedelame nic
  case 1: //computeResidueDissImplicit(w, g, BC, linSolver, setting);
    break;
  default: cout << "No such possibility for convection!" << endl;
    exit(15);
  }

}

#endif
