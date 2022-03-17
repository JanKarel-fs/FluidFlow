#ifndef VENKATAKRISHNAN_HPP
#define VENKATAKRISHNAN_HPP

#include <vector>
#include "cellField.hpp"
#include "grid.hpp"
#include "../geometry/vector.hpp"
#include <omp.h>

using namespace std;

template <typename var>
void venkatakrishnan(const CellField<var>& w, const CellField<Vector2<var> >& gradW,
		    const Grid& g, CellField<var>& psi) {

#pragma omp parallel for
  for (int i=psi.Imin(); i<psi.Imax(); i++) {
    for (int j=psi.Jmin(); j<psi.Jmax(); j++) {
      psi[i][j].zero();
    }
  }

#pragma omp parallel for
  for (int i=0; i<psi.M(); i++) {
    for (int j=0; j<psi.N(); j++) {
      psi[i][j].one();
    }
  }

#pragma omp parallel for
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

      double K = 5.;
      double dh = sqrt(g.volume(i, j));
      double eps2 = pow(K * dh, 3);
      
        
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
	      double delta1_max = w_max[m] - w[i][j][m]; 
	      double value = (((pow(delta1_max, 2) + eps2) * delta2
			       + 2.*pow(delta2, 2)*delta1_max) / (pow(delta1_max, 2)
			       + 2.*pow(delta2, 2) + delta1_max*delta2 + eps2)) / delta2;
	      psi_k = value;
	    }
	    if (delta2 < 0.) {
	      double delta1_min = w_min[m] - w[i][j][m];
	      double value = (((pow(delta1_min, 2) + eps2) * delta2
			       + 2.*pow(delta2, 2)*delta1_min) / (pow(delta1_min, 2)
			       + 2.*pow(delta2, 2) + delta1_min*delta2 + eps2)) / delta2;
	      psi_k = value;
	    }
	  }
            
	  psi_ij[m] = min(psi_ij[m], psi_k);
	}

      }
    }
  }
}

#endif
