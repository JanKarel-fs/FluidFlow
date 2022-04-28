#ifndef STEPEXPLICIT_HPP
#define STEPEXPLICIT_HPP

#include "fvm/grid.hpp"
#include "fvm/cellField.hpp"
#include "fvm/computeResidue.hpp"
#include "sources/typedefs.hpp"
#include "sources/settings.hpp"
#include "sources/setBoundaryConditions.hpp"
#include "sources/linearSolver.hpp"

template <typename var>
void stepExplicit(CellField<var>& w, CellField<var>& wOld, CellField<var>& res, const Grid& g,
		  const double& dt, const map<string, bcWithJacobian>& BC,
		  LinearSolver<var>& linSolver, const Settings& setting) {

  CellField<var> wStar(g);

  wStar = w;
  for (int k=0; k<setting.alphaRK.size(); k++) {
    setBoundaryConditions(wStar, g, setting, BC);

    computeResidue(wStar, g, res, setting);
      
    wStar = w + setting.alphaRK[k] * dt * res;
  }

  w = wStar;
  
}

#endif
