#ifndef SYMMETRYJACOBIAN_HPP
#define SYMMETRYJACOBIAN_HPP

#include <cmath>
#include "../geometry/vector.hpp"
#include "../geometry/matrix.hpp"
#include "../compressible.hpp"
#include "settings.hpp"

using namespace std;

Matrixd symmetryJacobian(const Compressible& wInside, const Vector2d& s, const Settings& setting);

#endif
