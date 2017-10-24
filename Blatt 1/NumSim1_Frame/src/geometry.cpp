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
	//multi_index_t *sizeptr = new multi_index_t(128, 128);
	//this-> _size = *sizeptr;
}


/// Loads a geometry from a file
void Geometry::Load(const char * file) {

}


/// Returns the number of cells in each dimension
const multi_index_t & Geometry::Size() const {
	// TODO: insert return statement here
}

/// Returns the length of the domain
const multi_real_t & Geometry::Length() const {
	// TODO: insert return statement here
}

/// Returns the meshwidth
const multi_real_t & Geometry::Mesh() const {
	// TODO: insert return statement here
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
