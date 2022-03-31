#ifndef GRID_XY_HPP
#define GRID_XY_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "grid.hpp"
#include "../geometry/vector.hpp"

using namespace std;

class Grid_xy : public Grid {
public:
  Grid_xy() {};
  Grid_xy(const string name1, const string name2, const int gh, const string& type);
  ~Grid_xy() {};

  void loadXYFiles(const string& name1, const string& name2,
		   vector<double>& xCoord, vector<double>& yCoord);
};

#endif
