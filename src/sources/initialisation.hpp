#ifndef INITIALISATION_HPP
#define INITIALISATION_HPP

#include "settings.hpp"
#include "../fvm/grid.hpp"
#include "../fvm/cellField.hpp"
#include "../compressible.hpp"
#include "../fvm/limiter.hpp"
#include "../fvm/zeroLimiter.hpp"
#include "../fvm/barthJespersen.hpp"
#include "../fvm/grad.hpp"
#include "../fvm/gradLSM.hpp"
#include "../fvm/zeroGrad.hpp"

void initialisation(CellField<Compressible>& w, const Grid& g, const Settings& setting);

#endif
