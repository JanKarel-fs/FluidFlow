#include "grid.hpp"

int Grid::Mnd() const {return Mnodes;}

int Grid::Nnd() const {return Nnodes;}

int Grid::Mvol() const {return Mvolumes;}

int Grid::Nvol() const {return Nvolumes;}

int Grid::ghostNum() const {return ghost;}

double Grid::x(const int i, const int j) const {
  return nodes[i][j].vertex.x;
}

double Grid::y(const int i, const int j) const {
  return nodes[i][j].vertex.y;
}

double Grid::volume(const int i, const int j) const {
  return volumes[i][j];
}

Point2d Grid::node(const int i, const int j) const {
  return nodes[i][j].vertex;
}

Point2d Grid::center(const int i, const int j) const {
  return centers[i][j];
}

Face Grid::faceI(const int i, const int j) const {
  return facesI[i][j];
}

Face Grid::faceJ(const int i, const int j) const {
  return facesJ[i][j];
}


void Face::updateI(const Grid& g, const int& i, const int& j) {
  const Point2d& a = g.node(i, j);
  const Point2d& b = g.node(i+1, j);

  center = 0.5 * (a + b);

  s = Vector2d(b.y-a.y, a.x-b.x);
}


void Face::updateJ(const Grid& g, const int& i, const int& j) {
  const Point2d& a = g.node(i, j);
  const Point2d& b = g.node(i, j+1);

  center = 0.5 * (a + b);

  s = Vector2d(b.y-a.y, a.x-b.x);
}


void Grid::update_faces() {
  for (int i=-ghost; i<Mvolumes+ghost; i++) {
    for (int j=-ghost; j<Nnodes+ghost; j++) {
      facesI[i][j].updateI(*this, i, j);
    }
  }

  for (int i=-ghost; i<Mnodes+ghost; i++) {
    for (int j=-ghost; j<Nvolumes+ghost; j++) {
      facesJ[i][j].updateJ(*this, i, j);
    }
  }
}
