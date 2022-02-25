#ifndef COMPUTERESIDUE_HPP
#define COMPUTERESIDUE_HPP

#include "grid.hpp"
#include "cellField.hpp"
#include "../geometry/vector.hpp"

template <typename var>
void computeResidue(const CellField<var>& w, const Grid& g, CellField<var>& res) {
  // vynulovani residua
  for (int i=0; i<res.M(); i++) {
    for (int j=0; j<res.N(); j++) {
      res[i][j].zero();
    }
  }
  
  // vypocet toku stenami ve smeru i
  for (int i=0; i<w.M(); i++) {
    for (int j=0; j<w.N()+1; j++) {
      
      var wl = w[i][j];
      var wr = w[i][j-1];
      // Vector2d f(g.node(i, j), g.node(i+1, j));
      // Vector2d s(f.y, -f.x);
      Vector2d s = g.faceI(i, j).s;

      var flx = var::flux(wl, wr, s);

      res[i][j] -= flx;
      res[i][j-1] += flx;
    }
  }

  // vypocet toku stenami ve smeru j
  for (int i=0; i<w.M()+1; i++) {
    for (int j=0; j<w.N(); j++) {
      
      var wl = w[i-1][j];
      var wr = w[i][j];
      // Vector2d f(g.node(i, j), g.node(i, j+1));
      // Vector2d s(f.y, -f.x);
      Vector2d s = g.faceJ(i, j).s;

      var flx = var::flux(wl, wr, s);

      res[i-1][j] -= flx;
      res[i][j] += flx;
    }
  }

  // deleni objemem bunky
  for (int i=0; i<res.M(); i++) {
    for (int j=0; j<res.N(); j++) {
      res[i][j] /= g.volume(i, j);
    }
  }

}

#endif
