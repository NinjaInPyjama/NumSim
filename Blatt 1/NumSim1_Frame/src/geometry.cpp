#include "geometry.hpp"

/// Constructs a default geometry:
// driven cavity with 128 x 128 grid, no-slip boundary conditions
// as shown below
//
//      u=1, v=0
//    -------------
//    |           |
// u=0|           |u=0
// v=0|           |v=0
//    |           |
//    |           |
//    -------------
//      u=0, v=0
Geometry::Geometry() {
	_size = multi_index_t(128); 
	_length = multi_real_t(1.0);
	_h = multi_real_t(1.0/128);

	_velocity = multi_real_t(1.0, 0.0);
	_pressure = 0.0;
}


/// Loads a geometry from a file
void Geometry::Load(const char * file) {
	// TODO
}


/// Returns the number of cells in each dimension
const multi_index_t & Geometry::Size() const {
	return _size;
}

/// Returns the length of the domain
const multi_real_t & Geometry::Length() const {
	return _length;
}

/// Returns the meshwidth
const multi_real_t & Geometry::Mesh() const {
	return _h;
}


/// Updates the velocity field u
void Geometry::Update_U(Grid * u) const {
	BoundaryIterator bit = BoundaryIterator(new Geometry());
	bit.setBoundary(0);
	for(bit.First(); bit.Valid(); bit.Next()) {
		u->Cell(bit) = _velocity[0];
	}
	for(index_t i=1; i<4; i++) {
		bit.setBoundary(i);
		for(bit.First(); bit.Valid(); bit.Next()) {
			u->Cell(bit) = 0.0;
		}
	}
}

/// Updates the velocity field v
void Geometry::Update_V(Grid * v) const {
	BoundaryIterator bit = BoundaryIterator(new Geometry());
	bit.setBoundary(0);
	for(bit.First(); bit.Valid(); bit.Next()) {
		u->Cell(bit) = _velocity[1];
	}
	for(index_t i=1; i<4; i++) {
		bit.setBoundary(i);
		for(bit.First(); bit.Valid(); bit.Next()) {
			u->Cell(bit) = 0.0;
		}
	}
}

/// Updates the pressure field p
void Geometry::Update_P(Grid * p) const {

}