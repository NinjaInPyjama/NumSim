#include "iterator.hpp"

/// Default Constructor
Iterator::Iterator() {}

/// Constructs a new Iterator depending on a geometry
Iterator::Iterator(const Geometry *geom){
	this -> _geom = geom;
	this -> _value = 0;
	this -> _valid = true;
}

/// Constructs a new Iterator on a geometry with a defined starting value
Iterator::Iterator(const Geometry * geom, const index_t &value) {
	this -> _geom = geom;
	this -> _value = value;
	this -> _valid = true;
}


/// Returns the current position value
const index_t &Iterator::Value() const {
	return _value;
}

/// Cast operator to convert Iterators to integers
Iterator::operator const index_t&() const {
	return _value;
}

/// Returns the position coordinates
multi_index_t Iterator::Pos() const {
	if(_valid) return multi_index_t(0, 0);
	else return multi_index_t(-1);
}


/// Sets the iterator to the first element
void Iterator::First() {

}

/// Goes to the next element of the iterator, disables it if position is end
void Iterator::Next() {

}


/// Checks if the iterator still has a valid value
bool Iterator::Valid() const {
	return false;
}


/// Returns an Iterator that is located left from this one.
// if we are at the left boundary, the cell sees itself
Iterator Iterator::Left() const {
	return nullptr;
}

/// Returns an Iterator that is located right from this one
// If we are at the right boundary, the cell sees itself
Iterator Iterator::Right() const {
	return nullptr;
}

/// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
Iterator Iterator::Top() const {
	return nullptr;
}

/// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
Iterator Iterator::Down() const {
	return nullptr;
}


//------------------------------------------------------------------------------
/** Iterator for interior cells
*/
/// Construct a new InteriorIterator
InteriorIterator::InteriorIterator(const Geometry * geom) {

}


/// Sets the iterator to the first element
void InteriorIterator::First() {

}

/// Goes to the next element of the iterator, disables it if position is end
void InteriorIterator::Next() {

}


//------------------------------------------------------------------------------
/** Iterator for domain boundary cells.
*/
/// Constructs a new BoundaryIterator
BoundaryIterator::BoundaryIterator(const Geometry * geom) {

}


/// Sets the boundary to iterate
void BoundaryIterator::SetBoundary(const index_t & boundary) {

}


/// Sets the iterator to the first element
void BoundaryIterator::First() {

}

/// Goes to the next element of the iterator, disables it if position is end
void BoundaryIterator::Next() {

}
