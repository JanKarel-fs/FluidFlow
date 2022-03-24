#ifndef SETNODEFIELD_HPP
#define SETNODEFIELD_HPP

#include <vector>
#include "fvm/cellField.hpp"
#include "fvm/nodeField.hpp"
#include "fvm/grid.hpp"
#include <omp.h>

using namespace std;

template <typename var>
void setNodeField(const CellField<var>& w, const Grid& g, NodeField<var>& wNode) {

#pragma omp parallel for
  for (int i=0; i<wNode.M(); i++) {
    for (int j=0; j<wNode.N(); j++) {
      vector<var> wk(4);

      wk[0] = w[i-1][j];
      wk[1] = w[i][j];
      wk[2] = w[i-1][j-1];
      wk[3] = w[i][j-1];

      vector<double> alpha = g.alpha(i, j);

      wNode[i][j].zero();

      for (int k=0; k<wk.size(); k++) {
	wNode[i][j] += alpha[k] * wk[k];
      }

    }
  }
}

#endif
