#include "grid.hpp"

/// Constructs a grid based on a geometry
Grid::Grid(const Geometry * geom) {

}

/// Constructs a grid based on a geometry with an offset
// @param geom   Geometry information
// @param offset distance of staggered grid point to cell's anchor point;
//               (anchor point = lower left corner)
Grid::Grid(const Geometry * geom, const multi_real_t & offset) {

}


/// Deletes the grid
Grid::~Grid() {
	
}


///     Initializes the grid with a value
void Grid::Initialize(const real_t & value) {

}


/// Write access to the grid cell at position [it]
real_t & Grid::Cell(const Iterator & it) {
	// TODO: insert return statement here
}

/// Read access to the grid cell at position [it]
const real_t & Grid::Cell(const Iterator & it) const {
	// TODO: insert return statement here
}


/// Interpolate the value at a arbitrary position
real_t Grid::Interpolate(const multi_real_t & pos) const {
	return 0;
}


/// Computes the left-sided difference quatient in x-dim at [it]
real_t Grid::dx_l(const Iterator & it) const {
	return 0;
}

/// Computes the right-sided difference quatient in x-dim at [it]
real_t Grid::dx_r(const Iterator & it) const {
	return 0;
}

/// Computes the left-sided difference quatient in y-dim at [it]
real_t Grid::dy_l(const Iterator & it) const {
	return 0;
}

/// Computes the right-sided difference quatient in x-dim at [it]
real_t Grid::dy_r(const Iterator & it) const {
	return 0;
}

/// Computes the central difference quatient of 2nd order in x-dim at [it]
real_t Grid::dxx(const Iterator & it) const {
	return 0;
}

/// Computes the central difference quatient of 2nd order in y-dim at [it]
real_t Grid::dyy(const Iterator & it) const {
	return 0;
}


/// Computes u*du/dx with the donor cell method
real_t Grid::DC_udu_x(const Iterator & it, const real_t & alpha) const {
	return 0;
}

/// Computes v*du/dy with the donor cell method
real_t Grid::DC_vdu_y(const Iterator & it, const real_t & alpha, const Grid * v) const {
	return 0;
}

/// Computes u*dv/dx with the donor cell method
real_t Grid::DC_udv_x(const Iterator & it, const real_t & alpha, const Grid * u) const {
	return 0;
}

/// Computes v*dv/dy with the donor cell method
real_t Grid::DC_vdv_y(const Iterator & it, const real_t & alpha) const {
	return 0;
}


/// Returns the maximal value of the grid
real_t Grid::Max() const {
	return 0;
}

/// Returns the minimal value of the grid
real_t Grid::Min() const {
	return 0;
}

/// Returns the absolute maximal value
real_t Grid::AbsMax() const {
	return 0;
}

/// Returns a pointer to the raw data
real_t * Grid::Data() {
	return nullptr;
}
