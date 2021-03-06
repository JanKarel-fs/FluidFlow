#ifndef COMPUTERESIDUECONV_HPP
#define COMPUTERESIDUECONV_HPP

#include "grid.hpp"
#include "cellField.hpp"
#include "grad.hpp"
#include "limiter.hpp"
#include "../geometry/vector.hpp"
#include <omp.h>

template <typename var>
void computeResidueConv(const CellField<var>& w, const Grid& g, CellField<var>& res) {

  CellField<Vector2<var> > gradW(g);
  CellField<var> psi(g);
  grad<var>(w, g, gradW);
  limiter<var>(w, gradW, g, psi); 
  
  // vypocet toku stenami ve smeru i
#pragma omp parallel for
  for (int i=0; i<w.M(); i++) {
    for (int j=0; j<w.N()+1; j++) {
      const Face& f = g.faceI(i, j);

      Vector2d rL(g.center(i, j), f.center);
      Vector2d rR(g.center(i, j-1), f.center);
      
      var wl = w[i][j] + psi[i][j] * (gradW[i][j].x * rL.x + gradW[i][j].y * rL.y);
      var wr = w[i][j-1] + psi[i][j-1] * (gradW[i][j-1].x * rR.x + gradW[i][j-1].y * rR.y);

      var flx = var::flux(wl, wr, f.s);

      res[i][j] -= flx;
      res[i][j-1] += flx;
    }
  }

  // vypocet toku stenami ve smeru j
#pragma omp parallel for
  for (int j=0; j<w.N(); j++) {
    for (int i=0; i<w.M()+1; i++) {
      const Face& f = g.faceJ(i, j);

      Vector2d rL(g.center(i-1, j), f.center);
      Vector2d rR(g.center(i, j), f.center);
      
      var wl = w[i-1][j] + psi[i-1][j] * (gradW[i-1][j].x * rL.x + gradW[i-1][j].y * rL.y);
      var wr = w[i][j] + psi[i][j] * (gradW[i][j].x * rR.x + gradW[i][j].y * rR.y);

      var flx = var::flux(wl, wr, f.s);

      res[i-1][j] -= flx;
      res[i][j] += flx;
    }
  }
}

#endif
