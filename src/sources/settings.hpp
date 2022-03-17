#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include <set>
#include "../geometry/vector.hpp"
#include "loadDataFile.hpp"
#include "findSection.hpp"

using namespace std;

class Settings {
public:
  int grid_type;
  int mCells;
  int nCells;
  int ghostCells;
  double pInit;
  double rhoInit;
  Vector2d uInit;
  double alpha;
  double Ma2is;
  int flux;
  double kappa;
  double p0;
  double rho0;
  double R;
  double Pr;
  int stop;
  double CFL;
  int temporalOrder;
  int spatialOrder;
  int limiter;
  vector<double> alphaRK;
  map<string, string> usedBC;

  Settings() {};
  Settings(const string& name);
  ~Settings() {};
};

#endif
