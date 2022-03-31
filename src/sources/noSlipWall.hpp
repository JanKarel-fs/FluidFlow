#ifndef NOSLIPWALL_HPP
#define NOSLIPWALL_HPP

#include <cmath>
#include "../geometry/vector.hpp"
#include "../compressible.hpp"
#include "settings.hpp"

using namespace std;

Compressible noSlipWall(const Compressible& wInside, const Vector2d& s, const Settings& setting);

#endif
