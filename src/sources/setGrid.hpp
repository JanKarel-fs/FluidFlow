#ifndef SETGRID_HPP
#define SETGRID_HPP

#include <iostream>
#include <cstdlib>
#include "../fvm/grid.hpp"
#include "../fvm/grid_gamm.hpp"
#include "settings.hpp"

using namespace std;

void setGrid(Grid& g, const Settings& setting);

#endif
