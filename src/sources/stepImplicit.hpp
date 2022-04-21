#ifndef STEPEIMPLICIT_HPP
#define STEPIMPLICIT_HPP

#include "fvm/grid.hpp"
#include "fvm/cellField.hpp"
#include "fvm/computeResidueImplicit.hpp"
#include "sources/typedefs.hpp"
#include "sources/settings.hpp"
#include "sources/setBoundaryConditions.hpp"
#include "sources/linearSolver.hpp"

template <typename var>
void stepImplicit(CellField<var>& w, CellField<var>& res, const Grid& g, const double& dt,
		  const map<string, bcWithJacobian>& BC, LinearSolver<var>& linSolver,
		  const Settings& setting) {

  linSolver.reset();

  setBoundaryConditions(w, g, setting, BC);

  computeResidueImplicit(w, g, BC, dt, linSolver, setting);

  linSolver.solve();

  linSolver.getResults(res);

  w = w + res;

  res = res / dt;
  
}

#endif
