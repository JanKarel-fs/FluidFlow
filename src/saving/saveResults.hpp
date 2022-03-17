#ifndef SAVERESULTS_HPP
#define SAVERESULTS_HPP

#include <fstream>
#include "../fvm/grid.hpp"
#include "../fvm/cellField.hpp"
#include "../fvm/nodeField.hpp"
#include "../compressible.hpp"
#include <omp.h>

using namespace std;

void saveResults(const CellField<Compressible>& w, const Grid& g);

#endif
