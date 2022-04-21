#ifndef TIMEINCREMENTFIRSTORDER_HPP
#define TIMEINCREMENTFIRSTORDER_HPP

#include "sources/linearSolver.hpp"

template <typename var>

void timeIncrementFirstOrder(LinearSolver<var>& linSolver, const double& dt) {
  
  for (int i=0; i<linSolver.nmb; i++) {
    MatSetValue(linSolver.A, i, i, 1./dt, ADD_VALUES);
  }
}

#endif
