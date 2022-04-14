#ifndef BCLIST_HPP
#define BCLIST_HPP

#include <map>
#include <string>
#include "typedefs.hpp"
#include "inlet.hpp"
#include "outlet.hpp"
#include "slipWall.hpp"
#include "noSlipWall.hpp"
#include "symmetry.hpp"
#include "homogeneousNeumann.hpp"
#include "inletJacobian.hpp"
#include "outletJacobian.hpp"
#include "slipWallJacobian.hpp"
#include "noSlipWallJacobian.hpp"
#include "symmetryJacobian.hpp"
#include "homogeneousNeumannJacobian.hpp"

using namespace std;

map<string, bcWithJacobian> bcList =
    {condition("inlet", bcWithJacobian(inlet, inletJacobian)),
     condition("outlet", bcWithJacobian(outlet, outletJacobian)),
     condition("slipWall", bcWithJacobian(slipWall, slipWallJacobian)),
     condition("noSlipWall", bcWithJacobian(noSlipWall, noSlipWallJacobian)),
     condition("symmetry", bcWithJacobian(symmetry, symmetryJacobian)),
     condition("homogeneousNeumann", bcWithJacobian(homogeneousNeumann,
						    homogeneousNeumannJacobian))};

#endif
