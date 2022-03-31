#ifndef HOMOGENEOUSNEUMANN_HPP
#define HOMOGENEOUSNEUMANN_HPP

#include <cmath>
#include "../geometry/vector.hpp"
#include "../compressible.hpp"
#include "settings.hpp"

using namespace std;

Compressible homogeneousNeumann(const Compressible& wInside, const Vector2d& s, const Settings& setting);

#endif
