#ifndef SYMMETRY_HPP
#define SYMMETRY_HPP

#include <cmath>
#include "../geometry/vector.hpp"
#include "../compressible.hpp"
#include "settings.hpp"

using namespace std;

Compressible symmetry(const Compressible& wInside, const Vector2d& s, const Settings& setting);

#endif
