#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <map>
#include <string>
#include <cmath>
#include "../geometry/point.hpp"
#include "../geometry/vector.hpp"
#include "../geometry/field.hpp"

using namespace std;

class Grid;

typedef void (*coefficients)(const Grid& g);

class Node {
public:
  Point2d vertex;
  vector<double> alpha;    // koeficienty (vahy) pro vypocet hodnoty neznamych ve vrcholu
};

class Face {
public:
  Point2d center;
  Vector2d s;          // normalovy vektor ke stene, prenasobeny velikosti steny
  string name;

  void updateI(const Grid& g, const int& i, const int& j);
  void updateJ(const Grid& g, const int& i, const int& j);
};

class Grid {
protected:
  Field2<Node> nodes;
  Field2<Face> facesI;
  Field2<Face> facesJ;
  Field2<Point2d> centers;
  Field2<double> volumes;

  int Mnodes;      //pocet uzlu ve smeru i
  int Nnodes;      //pocet uzlu ve smeru j
  int Mvolumes;    //pocet bunek ve smeru i
  int Nvolumes;    //pocet bunek ve smeru j

  int ghost;       //pocet vrstev pomocnych bunek

public:
  Grid() {};
  ~Grid() {};
  
  int Mnd() const;
  int Nnd() const;
  int Mvol() const;
  int Nvol() const;

  int ghostNum() const;

  double x(const int i, const int j) const;
  double y(const int i, const int j) const;
  double volume(const int i, const int j) const;
  Point2d node(const int i, const int j) const;
  Point2d center(const int i, const int j) const;
  vector<double> alpha(const int i, const int j) const;
  Face faceI(const int i, const int j) const;
  Face faceJ(const int i, const int j) const;

  static void (*computeCoefficients)(const Grid& g);
  static void computeCoefficientsLSM(const Grid& g);
  static void computeCoefficientsWeight(const Grid& g);

  static map<string, coefficients> mCoefficients;

  void update_faces();
};

#endif
