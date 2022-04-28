#ifndef TIMEINCREMENTSECONDORDER_HPP
#define TIMEINCREMENTSECONDORDER_HPP

#include "sources/linearSolver.hpp"

template <typename var>
void timeIncrementSecondOrder(LinearSolver<var>& linSolver, const CellField<var>& w,
			      const CellField<var>& wOld, const double& dt) {
  
  for (int i=0; i<linSolver.nmb; i++) {
    MatSetValue(linSolver.A, i, i, 3./(2.*dt), ADD_VALUES);
  }

  for (int i=0; i<w.M(); i++) {
    for (int j=0; j<w.N(); j++) {
      int offset = var::nVars * (j*w.M() + i);
      for (int k=0; k<var::nVars; k++) {
	linSolver.rhs[offset + k] += 1./(2.*dt)*(w[i][j][k] - wOld[i][j][k]);
      }
    }
  }
}

#endif
