#ifndef COMPUTERESIDUECONVIMPLICIT_HPP
#define COMPUTERESIDUECONVIMPLICIT_HPP

#include "grid.hpp"
#include "cellField.hpp"
#include "grad.hpp"
#include "limiter.hpp"
#include "../geometry/vector.hpp"
#include "../geometry/matrix.hpp"
#include "sources/typedefs.hpp"
#include "sources/linearSolver.hpp"
#include <omp.h>

template <typename var>
void computeResidueConvImplicit(const CellField<var>& w, const Grid& g,
				const map<string, bcWithJacobian>& BC,
				LinearSolver<var>& linSolver, const Settings& setting) {

  CellField<Vector2<var> > gradW(g);
  CellField<var> psi(g);
  grad<var>(w, g, gradW);
  limiter<var>(w, gradW, g, psi); 
  
  // vypocet toku stenami ve smeru i
  //#pragma omp parallel for
  for (int i=0; i<w.M(); i++) {
    for (int j=0; j<w.N()+1; j++) {
      const Face& f = g.faceI(i, j);

      Vector2d rL(g.center(i, j), f.center);
      Vector2d rR(g.center(i, j-1), f.center);
      
      var wl = w[i][j] + psi[i][j] * (gradW[i][j].x * rL.x + gradW[i][j].y * rL.y);
      var wr = w[i][j-1] + psi[i][j-1] * (gradW[i][j-1].x * rR.x + gradW[i][j-1].y * rR.y);

      pair<pair<Matrixd, Matrixd>, var> increment = var::fluxImplicit(wl, wr, f.s);

      int leftOffset = var::nVars * (j*w.M() + i);
      int rightOffset = var::nVars * ((j-1)*w.M() + i);

      double leftVolume = g.volume(i, j);
      double rightVolume = g.volume(i, j-1);

      const var& flx = increment.second;
      pair<Matrixd, Matrixd>& Jacobians = increment.first;

      if (j==0 || j==w.N()) {
	var wb;
	int offset;
	if (j==0) {
	  wb = wl;
	  offset = leftOffset;
	}
	else {
	  wb = wr;
	  offset = rightOffset;
	}
	auto it = BC.find(f.name);
	Matrixd BJacobian = it->second.second(wb, f.s, setting);
	Jacobians.second = Jacobians.second * BJacobian;
	
	for (int k=0; k<var::nVars; k++) {
	  VecSetValue(linSolver.b, offset + k, -1.*flx[k] / leftVolume, ADD_VALUES);
	  for (int m=0; m<var::nVars; m++) {
	    MatSetValue(linSolver.A, offset + k, offset + m,
			Jacobians.first[k][m] / leftVolume, ADD_VALUES);
	    MatSetValue(linSolver.A, offset + k, offset + m,
			Jacobians.second[k][m] / leftVolume, ADD_VALUES);
	  }
	}
      }
      else {
	for (int k=0; k<var::nVars; k++) {
	  VecSetValue(linSolver.b, leftOffset + k, -1.*flx[k] / leftVolume, ADD_VALUES);
	  VecSetValue(linSolver.b, rightOffset + k, flx[k] / rightVolume, ADD_VALUES);
	  for (int m=0; m<var::nVars; m++) {
	    MatSetValue(linSolver.A, leftOffset + k, leftOffset + m,
			Jacobians.first[k][m] / leftVolume, ADD_VALUES);
	    MatSetValue(linSolver.A,leftOffset + k, rightOffset + m,
			Jacobians.second[k][m] / leftVolume, ADD_VALUES);
	
	    MatSetValue(linSolver.A, rightOffset + k, leftOffset + m,
			-1.*Jacobians.first[k][m] / rightVolume, ADD_VALUES);
	    MatSetValue(linSolver.A, rightOffset + k, rightOffset + m,
			-1.*Jacobians.second[k][m] / rightVolume, ADD_VALUES);
	  }
	}
      }
    }
  }

  // vypocet toku stenami ve smeru j
  //#pragma omp parallel for
  for (int j=0; j<w.N(); j++) {
    for (int i=0; i<w.M()+1; i++) {
      const Face& f = g.faceJ(i, j);

      Vector2d rL(g.center(i-1, j), f.center);
      Vector2d rR(g.center(i, j), f.center);
      
      var wl = w[i-1][j] + psi[i-1][j] * (gradW[i-1][j].x * rL.x + gradW[i-1][j].y * rL.y);
      var wr = w[i][j] + psi[i][j] * (gradW[i][j].x * rR.x + gradW[i][j].y * rR.y);
      
      pair<pair<Matrixd, Matrixd>, var> increment = var::fluxImplicit(wl, wr, f.s);

      int leftOffset = var::nVars * (j*w.M() + i-1);
      int rightOffset = var::nVars * (j*w.M() + i);

      double leftVolume = g.volume(i-1, j);
      double rightVolume = g.volume(i, j);

      const var& flx = increment.second;
      pair<Matrixd, Matrixd>& Jacobians = increment.first;

      if (i==0 || i==w.M()) {
	var wb;
	int offset;
	if (i==0) {
	  wb = wr;
	  offset = rightOffset;
	}
	else {
	  wb = wl;
	  offset = leftOffset;
	}
	auto it = BC.find(f.name);
	Matrixd BJacobian = it->second.second(wb, f.s, setting);
	Jacobians.second = Jacobians.second * BJacobian;
	
	for (int k=0; k<var::nVars; k++) {
	  VecSetValue(linSolver.b, offset + k, -1.*flx[k] / leftVolume, ADD_VALUES);
	  for (int m=0; m<var::nVars; m++) {
	    MatSetValue(linSolver.A, offset + k, offset + m,
			Jacobians.first[k][m] / leftVolume, ADD_VALUES);
	    MatSetValue(linSolver.A, offset + k, offset + m,
			Jacobians.second[k][m] / leftVolume, ADD_VALUES);
	  }
	}
      }
      else {
	for (int k=0; k<var::nVars; k++) {
	  VecSetValue(linSolver.b, leftOffset + k, -1.*flx[k] / leftVolume, ADD_VALUES);
	  VecSetValue(linSolver.b, rightOffset + k, flx[k] / rightVolume, ADD_VALUES);
	  for (int m=0; m<var::nVars; m++) {
	    MatSetValue(linSolver.A, leftOffset + k, leftOffset + m,
			Jacobians.first[k][m] / leftVolume, ADD_VALUES);
	    MatSetValue(linSolver.A,leftOffset + k, rightOffset + m,
			Jacobians.second[k][m] / leftVolume, ADD_VALUES);
	
	    MatSetValue(linSolver.A, rightOffset + k, leftOffset + m,
			-1.*Jacobians.first[k][m] / rightVolume, ADD_VALUES);
	    MatSetValue(linSolver.A, rightOffset + k, rightOffset + m,
			-1.*Jacobians.second[k][m] / rightVolume, ADD_VALUES);
	  }
	}
      }
    }
  }
}

#endif
