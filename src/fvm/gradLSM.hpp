#ifndef GRADLSM_HPP
#define GRADLSM_HPP

#include <vector>
#include <cmath>
#include "grid.hpp"
#include "cellField.hpp"
#include "../geometry/vector.hpp"
#include "../geometry/point.hpp"
#include <omp.h>

template <typename var>
void gradLSM(const CellField<var>& w, const Grid& g, CellField<Vector2<var> >& gradW) {

#pragma omp parallel for
  for (int i=gradW.Imin(); i<gradW.Imax(); i++) {
    for (int j=gradW.Jmin(); j<gradW.Jmax(); j++) {
      gradW[i][j].x.zero();
      gradW[i][j].y.zero();
    }
  }

#pragma omp parallel for
  for (int i=0; i<g.Mvol(); i++) {
    for (int j=0; j<g.Nvol(); j++) {
      vector<Point2d> centers;
      vector<var> Neighbours;

      for (int p=i-1; p<=i+1; p++) {
	for (int r=j-1; r<=j+1; r++) {
	  if (p!=i || r!=j) {
	    centers.push_back(g.center(p, r));
	    Neighbours.push_back(w[p][r]);
	  }
	}
      }

      Point2d centerIJ = g.center(i, j);
          
      double I_xx = 0.;
      double I_yy = 0.;
      double I_xy = 0.;
      var J_x;
      J_x.zero();
      var J_y;
      J_y.zero();
         
      for (int k=0; k<centers.size(); k++) {
	double Jx = (centers[k].x - centerIJ.x);
	double Jy = (centers[k].y - centerIJ.y);
             
	I_xx += pow(Jx, 2);
	I_yy += pow(Jy, 2);
	I_xy += (Jx * Jy);
	J_x += (Jx * (Neighbours[k] - w[i][j]));
	J_y += (Jy * (Neighbours[k] - w[i][j]));
      }
         
      double D = I_xx * I_yy - pow(I_xy, 2);
         
      gradW[i][j].x = (I_yy * J_x - I_xy * J_y) / D;
      gradW[i][j].y = (I_xx * J_y - I_xy * J_x) / D;
    }
  }
}

#endif
