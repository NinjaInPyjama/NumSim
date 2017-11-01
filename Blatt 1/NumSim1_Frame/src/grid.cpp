#include "grid.hpp"
#include <iostream>
using namespace std;

/// Constructs a grid based on a geometry
Grid::Grid(const Geometry * geom) {
	_data = new real_t(geom->Size()[0] * geom->Size()[1]);
	_offset = 0.0;
	_geom = geom;
}

/// Constructs a grid based on a geometry with an offset
// @param geom   Geometry information
// @param offset distance of staggered grid point to cell's anchor point;
//               (anchor point = lower left corner)
Grid::Grid(const Geometry * geom, const multi_real_t & offset) {
	_data = new real_t(geom->Size()[0] * geom->Size()[1]);
	_offset = offset;
	_geom = geom;
}


/// Deletes the grid
Grid::~Grid() {}


///     Initializes the grid with a value
void Grid::Initialize(const real_t & value) {
    index_t num_cells = _geom->Size()[0] * _geom->Size()[1];
    for (index_t i = 0; i < num_cells; i++) {
		_data[i] = value;
	}
}


/// Write access to the grid cell at position [it]
real_t & Grid::Cell(const Iterator & it) {
	return _data[it.Value()];
}

/// Read access to the grid cell at position [it]
const real_t & Grid::Cell(const Iterator & it) const {
	return _data[it.Value()];
}


/// Interpolate the value at a arbitrary position
real_t Grid::Interpolate(const multi_real_t & pos) const {
    index_t pos_x = index_t(pos[0] - _offset[0]);
    //cout << pos_x << endl;
    index_t pos_y = index_t(pos[1] - _offset[1]);
    //cout << pos_y << endl;
    
    real_t unten_links = _data[pos_x + pos_y*_geom->Size()[0]];
    //cout << unten_links << endl;
    real_t unten_rechts = _data[pos_x + 1 + pos_y*_geom->Size()[0]];
    //cout << unten_rechts << endl;
    real_t oben_links = _data[pos_x + (pos_y + 1)*_geom->Size()[0]];
    //cout << oben_links << endl;
    real_t oben_rechts = _data[pos_x + 1 + (pos_y + 1)*_geom->Size()[0]];
    //cout << oben_rechts << endl;
    
    real_t anteil_x = pos[0] - _offset[0] - real_t(pos_x);
    //cout << anteil_x << endl;
    real_t anteil_y = pos[1] - _offset[1] - real_t(pos_y);
    //cout << anteil_y << endl;
	return (unten_links*(1.0 - anteil_x) + anteil_x*unten_rechts)*(1.0-anteil_y) + anteil_y*( oben_links*(1.0 - anteil_x) + anteil_x*oben_rechts );
}


/// Computes the left-sided difference quatient in x-dim at [it]
real_t Grid::dx_l(const Iterator & it) const {
	return (Cell(it.Left()) - Cell(it))/_geom->Mesh()[0];
}

/// Computes the right-sided difference quatient in x-dim at [it]
real_t Grid::dx_r(const Iterator & it) const {
	return (Cell(it.Right()) - Cell(it)) / _geom->Mesh()[0];
}

/// Computes the left-sided difference quatient in y-dim at [it]
real_t Grid::dy_l(const Iterator & it) const {
	return (Cell(it.Top()) - Cell(it)) / _geom->Mesh()[1];
}

/// Computes the right-sided difference quatient in x-dim at [it]
real_t Grid::dy_r(const Iterator & it) const {
	return (Cell(it.Down()) - Cell(it)) / _geom->Mesh()[1];
}

/// Computes the central difference quatient of 2nd order in x-dim at [it]
real_t Grid::dxx(const Iterator & it) const {
	return (Cell(it.Right()) - Cell(it.Left())) / (2 * _geom->Mesh()[0]);
}

/// Computes the central difference quatient of 2nd order in y-dim at [it]
real_t Grid::dyy(const Iterator & it) const {
	return (Cell(it.Down()) - Cell(it.Top())) / (2 * _geom->Mesh()[1]);
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
	real_t max = _data[0];
	index_t num_cells = _geom->Size()[0] * _geom->Size()[1];
	for (index_t i = 1; i < num_cells; i++) {
		if(max < _data[i]) max = _data[i];
	}
	return max;
}

/// Returns the minimal value of the grid
real_t Grid::Min() const {
	real_t min = _data[0];
	index_t num_cells = _geom->Size()[0] * _geom->Size()[1];
	for (index_t i = 1; i < num_cells; i++) {
		if (min > _data[i]) min = _data[i];
	}
	return min;
}

/// Returns the absolute maximal value
real_t Grid::AbsMax() const {
	real_t max = abs(_data[0]);
	index_t num_cells = _geom->Size()[0] * _geom->Size()[1];
	for (index_t i = 1; i < num_cells; i++) {
		if (max < abs(_data[i])) max = abs(_data[i]);
	}
	return max;
}

/// Returns a pointer to the raw data
real_t * Grid::Data() {
	return _data;
}
