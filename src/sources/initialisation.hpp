#ifndef INITIALISATION_HPP
#define INITIALISATION_HPP

#include "settings.hpp"
#include "../fvm/grid.hpp"
#include "../fvm/cellField.hpp"
#include "../compressible.hpp"

void initialisation(CellField<Compressible>& w, const Grid& g, const Settings& setting);

#endif
