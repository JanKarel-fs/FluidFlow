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

vector<double> Grid::alpha(const int i, const int j) const {
  return nodes[i][j].alpha;
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


void Grid::computeCoefficientsWeight(const Grid& g) {
  for (int i=0; i<g.Mnd(); i++) {
    for (int j=0; j<g.Nnd(); j++) {
      double volume = g.volume(i-1, j) + g.volume(i, j) + g.volume(i-1, j-1) + g.volume(i, j-1);

      g.nodes[i][j].alpha.resize(4);
      g.nodes[i][j].alpha[0] = (g.volume(i-1, j) / volume);
      g.nodes[i][j].alpha[1] = (g.volume(i, j) / volume);
      g.nodes[i][j].alpha[2] = (g.volume(i-1, j-1) / volume);
      g.nodes[i][j].alpha[3] = (g.volume(i, j-1) / volume);
    }
  }
}

void Grid::computeCoefficientsLSM(const Grid& g) {
  for (int i=0; i<g.Mnd(); i++) {
    for (int j=0; j<g.Nnd(); j++) {
      double Ixx = 0.;
      double Iyy = 0.;
      double Ixy = 0.;
      double Rx = 0.;
      double Ry = 0;

      vector<Point2d> neighbours(4);
      neighbours[0] = g.center(i-1, j);
      neighbours[1] = g.center(i, j);
      neighbours[2] = g.center(i-1, j-1);
      neighbours[3] = g.center(i, j-1);

      Point2d V = g.node(i, j);

      for (int k=0; k<neighbours.size(); k++) {
	double Jx = neighbours[k].x - V.x;
	double Jy = neighbours[k].y - V.y;

	Ixx += Jx*Jx;
	Iyy += Jy*Jy;
	Ixy += Jx*Jy;
	Rx += Jx;
	Ry += Jy;
      }

      double D = Ixx*Iyy - pow(Ixy, 2);

      double lx = (Ry*Ixy - Rx*Iyy) / D;
      double ly = (Rx*Ixy - Ry*Ixx) / D;
      
      int N = neighbours.size();
      g.nodes[i][j].alpha.resize(N);

      for (int k=0; k<neighbours.size(); k++) {
	g.nodes[i][j].alpha[k] = (1. + lx * (neighbours[k].x - V.x) + ly * (neighbours[k].y - V.y)) / (N + lx*Rx + ly*Ry);
      }
      
    }
  }
}

map<string, coefficients> Grid::mCoefficients = {pair<string, coefficients>("Weight", Grid::computeCoefficientsWeight),
                                                 pair<string, coefficients>("LSM", Grid::computeCoefficientsLSM)};

void (*Grid::computeCoefficients)(const Grid& g);

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
