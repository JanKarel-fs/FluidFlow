#ifndef OUTLETJACOBIAN_HPP
#define OUTLETJACOBIAN_HPP

#include <cmath>
#include "../geometry/vector.hpp"
#include "../geometry/matrix.hpp"
#include "../compressible.hpp"
#include "settings.hpp"

using namespace std;

Matrixd outletJacobian(const Compressible& wInside, const Vector2d& s, const Settings& setting);

#endif
