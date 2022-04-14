#ifndef TYPEDEFS_HPP
#define TYPEDEFS_HPP

#include <utility>
#include "../compressible.hpp"
#include "../geometry/vector.hpp"
#include "../geometry/matrix.hpp"
#include "settings.hpp"

using namespace std;

typedef Compressible (*bCondition)(const Compressible& wInside, const Vector2d& s, const Settings& setting);

typedef Matrixd (*bJacobian)(const Compressible& wInside, const Vector2d& s, const Settings& setting);

typedef pair<bCondition, bJacobian> bcWithJacobian;

typedef pair<string, bcWithJacobian> condition;

#endif
