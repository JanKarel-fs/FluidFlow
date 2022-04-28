#ifndef STEP_HPP
#define STEP_HPP

#include "fvm/grid.hpp"
#include "fvm/cellField.hpp"
#include "sources/typedefs.hpp"
#include "sources/settings.hpp"
#include "sources/linearSolver.hpp"

using namespace std;

template <typename var>
void (*step)(CellField<var>& w, CellField<var>& wOld, CellField<var>& res, const Grid& g,
	     const double& dt, const map<string, bcWithJacobian>& BC, LinearSolver<var>& linSolver,
	     const Settings& setting);

#endif
