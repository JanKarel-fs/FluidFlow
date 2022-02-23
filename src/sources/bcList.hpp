#ifndef BCLIST_HPP
#define BCLIST_HPP

#include <map>
#include <string>
#include "typedefs.hpp"
#include "inlet.hpp"
#include "outlet.hpp"
#include "slipWall.hpp"

using namespace std;

map<string, bCondition> bcList = {condition("inlet", inlet),
				  condition("outlet", outlet),
				  condition("slipWall", slipWall)};

#endif
