#ifndef BARTHJESPERSEN_HPP
#define BARTHJESPERSEN_HPP

#include <vector>
#include "cellField.hpp"
#include "grid.hpp"
#include "../geometry/vector.hpp"

using namespace std;

template <typename var>
void barthJespersen(const CellField<var>& w, const CellField<Vector2<var> >& gradW,
		    const Grid& g, CellField<var>& psi) {

  for (int i=psi.Imin(); i<psi.Imax(); i++) {
    for (int j=psi.Jmin(); j<psi.Jmax(); j++) {
      psi[i][j].zero();
    }
  }

  for (int i=0; i<psi.M(); i++) {
    for (int j=0; j<psi.N(); j++) {
      psi[i][j].one();
    }
  }

  for (int i=0; i<w.M(); i++) {
    for (int j=0; j<w.N(); j++) {
      Point2d centerIJ = g.center(i, j);

      vector<var> Neighbours(4);
      vector<Point2d> centers(4);

      Neighbours[0] = w[i-1][j]; Neighbours[1] = w[i+1][j] ;
      Neighbours[2] = w[i][j-1]; Neighbours[3] = w[i][j+1];

      centers[0] = g.faceI(i, j).center;
      centers[1] = g.faceJ(i+1, j).center;
      centers[2] = g.faceI(i, j+1).center;
      centers[3] = g.faceJ(i, j).center;
 
      var& psi_ij = psi[i][j];
      var w_max, w_min;
      w_max = w_min = w[i][j];
        
      for (int k=0; k<Neighbours.size(); k++) {
	w_max = var::max(w_max, Neighbours[k]);
	w_min = var::min(w_min, Neighbours[k]);
      }
      
        
      for (int k=0; k<centers.size(); k++) {   
	double delta2;
	double psi_k;
	Vector2d r(g.center(i, j), centers[k]);

	for (int m=0; m<var::nVars; m++) {
	  delta2 = (gradW[i][j].x[m] * r.x + gradW[i][j].y[m] * r.y);
            
	  if (fabs(delta2) < 1e-8) {
	    psi_k = 1.;
	  }
	  else {
	    if (delta2 > 0.) {
	      double value = (w_max[m] - w[i][j][m]) / delta2;
	      psi_k = min(1., value);
	    }
	    if (delta2 < 0.) {
	      double value = (w_min[m] - w[i][j][m]) / delta2;
	      psi_k = min(1., value);
	    }
	  }
            
	  psi_ij[m] = min(psi_ij[m], psi_k);
	}

      }
    }
  }
}

#endif
