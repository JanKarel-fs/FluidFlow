#ifndef STORENORMRESIDUE_HPP
#define STORENORMRESIDUE_HPP

#include <cmath>
#include <fstream>
#include "../fvm/grid.hpp"
#include "../fvm/cellField.hpp"
#include "../compressible.hpp"

using namespace std;

void storeNormResidue(const CellField<Compressible>& res, const Grid& g, const int& it);

#endif
