#ifndef ZEROGRAD_HPP
#define ZEROGRAD_HPP

#include "grid.hpp"
#include "cellField.hpp"
#include <omp.h>

template <typename var>
void zeroGrad(const CellField<var>& w, const Grid& g, CellField<Vector2<var> >& gradW) {

#pragma omp parallel for
  for (int i=gradW.Imin(); i<gradW.Imax(); i++) {
    for (int j=gradW.Jmin(); j<gradW.Jmax(); j++) {
      gradW[i][j].x.zero();
      gradW[i][j].y.zero();
    }
  }
  
}

#endif
