#include "grid_gamm.hpp"

Grid_gamm::Grid_gamm(const int m, const int n, const int gh) {
  //m, n pocet bunek ve smeru i, j
  Mnodes = m+1;
  Nnodes = n+1;
  Mvolumes = m;
  Nvolumes = n;
  ghost = gh;

  nodes.allocate(-ghost, Mnodes+ghost, -ghost, Nnodes+ghost);
  centers.allocate(-ghost, Mvolumes+ghost, -ghost, Nvolumes+ghost);
  volumes.allocate(-ghost, Mvolumes+ghost, -ghost, Nvolumes+ghost);
  facesI.allocate(-ghost, Mvolumes+ghost, -ghost, Nnodes+ghost);
  facesJ.allocate(-ghost, Mnodes+ghost, -ghost, Nvolumes+ghost);

  //generovani souradnic

  double dx = 3. / Mvolumes;
  
  for (int i=0; i<Mnodes; i++) {
    double y0 = 0.;
    double y1 = 1.;
    double x = i*dx;
    if (x>1. && x<2.) y0 = sqrt(1.3*1.3 - (x-1.5)*(x-1.5)) - 1.2;
    double dy = (y1 - y0) / Nvolumes;
    for (int j=0; j<Nnodes; j++) {
      double y = y0 + j*dy;
      nodes[i][j].vertex = Point2d(x, y);
    }
  }

  // vytvoreni pomocnych bunek
  for (int j=0; j<Nnodes; j++) {
    for (int k=1; k<=ghost; k++) {
      //leva hranice
      Vector2d s(nodes[1][j].vertex, nodes[0][j].vertex);
       nodes[-k][j].vertex = nodes[0][j].vertex + k * s;

      //prava hranice
      s = Vector2d(nodes[Mnodes-2][j].vertex, nodes[Mnodes-1][j].vertex);
      nodes[Mnodes-1+k][j].vertex = nodes[Mnodes-1][j].vertex + k * s;
    }
  }

  for (int i=-ghost; i<Mnodes+ghost; i++) {
    for (int k=1; k<=ghost; k++) {
      //dolni hranice
      Vector2d s(nodes[i][1].vertex, nodes[i][0].vertex);
      nodes[i][-k].vertex = nodes[i][0].vertex + k * s;

      //horni hranice
      s = Vector2d(nodes[i][Nnodes-2].vertex, nodes[i][Nnodes-1].vertex);
      nodes[i][Nnodes-1+k].vertex = nodes[i][Nnodes-1].vertex + k * s;
    }
  }

  //vypocet objemu bunek a stredu bunek
  for (int i=-ghost; i<Mvolumes+ghost; i++) {
    for (int j=-ghost; j<Nvolumes+ghost; j++) {
      Vector2d u1(nodes[i][j].vertex, nodes[i+1][j+1].vertex);
      Vector2d u2(nodes[i+1][j].vertex, nodes[i][j+1].vertex);

      volumes[i][j] = 0.5 * fabs(cross(u1, u2));
      centers[i][j] = 0.25 * (node(i, j) + node(i+1, j) + node(i+1, j+1) + node(i, j+1));
    }
  }

  update_faces();

  for (int i=-ghost; i<Mvolumes+ghost; i++) {
    for (int j=-ghost; j<Nnodes+ghost; j++) {
      facesI[i][j].name = "internal";

      if (j == 0 || j == (Nnodes-1)) facesI[i][j].name = "wall";
    }
  }

  for (int i=-ghost; i<Mnodes+ghost; i++) {
    for (int j=-ghost; j<Nvolumes+ghost; j++) {
      facesJ[i][j].name = "internal";

      if (i == 0) facesJ[i][j].name = "inlet";
      if (i == (Mnodes-1)) facesJ[i][j].name = "outlet";
    }
  }
}
