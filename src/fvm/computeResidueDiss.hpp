#ifndef COMPUTERESIDUEDISS_HPP
#define COMPUTERESIDUEDISS_HPP

#include "grid.hpp"
#include "cellField.hpp"
#include "nodeField.hpp"
#include "grad.hpp"
#include "limiter.hpp"
#include "../geometry/vector.hpp"
#include "sources/setNodeField.hpp"
#include "firstDerivatives.hpp"
#include "primitiveVars.hpp"
#include <omp.h>

template <typename var>
void computeResidueDiss(const CellField<var>& w, const Grid& g, CellField<var>& res) {

  CellField<PrimitiveVars> pVars;
  NodeField<PrimitiveVars> pVarsNode;

#pragma omp parallel for
  for (int i=w.Imin(); i<w.Imax(); i++) {
    for (int j=w.Jmin(); j<w.Jmax(); j++) {
      pVars[i][j] = PrimitiveVars::set(w[i][j]);
    }
  }

  setNodeField(pVars, g, pVarsNode);
  
  // vypocet toku stenami ve smeru i
#pragma omp parallel for
  for (int i=0; i<w.M(); i++) {
    for (int j=0; j<w.N()+1; j++) {
      const Face& f = g.faceI(i, j);

      const PrimitiveVars& pVarsL = pVars[i][j];
      const PrimitiveVars& pVarsR = pVars[i][j-1];
      const PrimitiveVars& pVarsA = pVarsNode[i][j];
      const PrimitiveVars& pVarsB = pVarsNode[i+1][j];

      Point2d L = g.center(i, j);
      Point2d R = g.center(i, j-1);
      Point2d A = g.node(i, j);
      Point2d B = g.node(i+1, j);

      Vector2<PrimitiveVars> grad_pVars;
      firstDerivatives(pVarsL, pVarsR, pVarsA, pVarsB, L, R, A, B, f.s, grad_pVars);

      PrimitiveVars pVarsFace = 0.5 * (pVarsL + pVarsR);
      var wFace = 0.5 * (w[i][j] + w[i][j-1]);

      var flx = var::fluxDissipative(wFace, pVarsFace, grad_pVars, f.s);

      res[i][j] += flx;
      res[i][j-1] -= flx;
    }
  }

  // vypocet toku stenami ve smeru j
#pragma omp parallel for
  for (int j=0; j<w.N(); j++) {
    for (int i=0; i<w.M()+1; i++) {
      const Face& f = g.faceJ(i, j);

      const PrimitiveVars& pVarsL = pVars[i-1][j];
      const PrimitiveVars& pVarsR = pVars[i][j];
      const PrimitiveVars& pVarsA = pVarsNode[i][j];
      const PrimitiveVars& pVarsB = pVarsNode[i][j+1];

      Point2d L = g.center(i-1, j);
      Point2d R = g.center(i, j);
      Point2d A = g.node(i, j);
      Point2d B = g.node(i, j+1);

      Vector2<PrimitiveVars> grad_pVars;
      firstDerivatives(pVarsL, pVarsR, pVarsA, pVarsB, L, R, A, B, f.s, grad_pVars);

      PrimitiveVars pVarsFace = 0.5 * (pVarsL + pVarsR);
      var wFace = 0.5 * (w[i-1][j] + w[i][j]);

      var flx = var::fluxDissipative(wFace, pVarsFace, grad_pVars, f.s);

      res[i-1][j] += flx;
      res[i][j] -= flx;
    }
  }
}

#endif
