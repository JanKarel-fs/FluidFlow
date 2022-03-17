#ifndef GRID_GAMM_HPP
#define GRID_GAMM_HPP

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "grid.hpp"
#include "../geometry/vector.hpp"

class Grid_gamm : public Grid {
public:
  Grid_gamm() {};
  Grid_gamm(const int m, const int n, const int gh, const string& type);
  ~Grid_gamm() {};
};

#endif
