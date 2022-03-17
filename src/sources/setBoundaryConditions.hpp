#ifndef SETBOUNDARYCONDITIONS_HPP
#define SETBOUNDARYCONDITIONS_HPP

#include <map>
#include <string>
#include "settings.hpp"
#include "../geometry/vector.hpp"
#include "../fvm/grid.hpp"
#include "../fvm/cellField.hpp"
#include "../compressible.hpp"
#include "typedefs.hpp"
#include <omp.h>

using namespace std;

void setBoundaryConditions(CellField<Compressible>& w, const Grid& g, const Settings& setting,
			   const map<string, bCondition>& BC);

#endif
