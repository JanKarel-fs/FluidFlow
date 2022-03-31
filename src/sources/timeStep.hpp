#ifndef TIMESTEP_HPP
#define TIMESTEP_HPP

#include <limits>
#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "../fvm/grid.hpp"
#include "../fvm/cellField.hpp"
#include "../geometry/vector.hpp"
#include "../compressible.hpp"
#include "settings.hpp"

double timeStep(const CellField<Compressible>& w, const Grid& g, const Settings& setting);

#endif
