#ifndef LIMITER_HPP
#define LIMITER_HPP
 
#include <vector>
#include "cellField.hpp"
#include "grid.hpp"

template <typename var>
void (*limiter)(const CellField<var>& w, const CellField<Vector2<var> >& gradW,
		const Grid& g, CellField<var>& psi);

#endif
