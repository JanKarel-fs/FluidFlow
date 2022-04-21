#ifndef TIMEINCREMENT_HPP
#define TIMEINCReMENT_HPP

#include "sources/linearSolver.hpp"

template <typename var>
void (*timeIncrement)(LinearSolver<var>& linSolver, const double& dt);

#endif
