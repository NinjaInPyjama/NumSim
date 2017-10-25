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
	this -> _size = multi_index_t(128); 
	this -> _length = multi_real_t(1.0);
	this -> _h = multi_real_t(1.0/128);

	this -> _velocity = multi_real_t(1.0, 0.0);
	this -> _pressure = 0.0;
}


/// Loads a geometry from a file
void Geometry::Load(const char * file) {
	// TODO
}


/// Returns the number of cells in each dimension
const multi_index_t & Geometry::Size() const {
	return this -> _size;
}

/// Returns the length of the domain
const multi_real_t & Geometry::Length() const {
	return this -> _length;
}

/// Returns the meshwidth
const multi_real_t & Geometry::Mesh() const {
	return this -> _h;
}


/// Updates the velocity field u
void Geometry::Update_U(Grid * u) const {

}

/// Updates the velocity field v
void Geometry::Update_V(Grid * v) const {

}

/// Updates the pressure field p
void Geometry::Update_P(Grid * p) const {

}
