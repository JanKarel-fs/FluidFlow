#ifndef GRAD_HPP
#define GRAD_HPP

#include "grid.hpp"
#include "cellField.hpp"

template <typename var>
void (*grad)(const CellField<var>& w, const Grid& g, CellField<Vector2<var> >& gradW);

#endif
