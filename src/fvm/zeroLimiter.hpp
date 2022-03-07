#ifndef ZEROLIMITER_HPP
#define ZEROLIMITER_HPP

#include "cellField.hpp"
#include "grid.hpp"

template <typename var>
void zeroLimiter(const CellField<var>& w, const CellField<Vector2<var> >& gradW,
		    const Grid& g, CellField<var>& psi) {

  for (int i=psi.Imin(); i<psi.Imax(); i++) {
    for (int j=psi.Jmin(); j<psi.Jmax(); j++) {
      psi[i][j].zero();
    }
  }
}

#endif
