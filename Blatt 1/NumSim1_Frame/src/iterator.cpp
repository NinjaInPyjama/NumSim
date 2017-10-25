#include "iterator.hpp"

/// Default Constructor
Iterator::Iterator() {}

/// Constructs a new Iterator depending on a geometry
Iterator::Iterator(const Geometry *geom){
	_geom = geom;
	First();
	_valid = true;
}

/// Constructs a new Iterator on a geometry with a defined starting value
Iterator::Iterator(const Geometry * geom, const index_t &value) {
	_geom = geom;
	_value = value;
	_valid = true;
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
	if (_valid) {
		index_t xPos = (_value % _geom->Size()[0]) + 1; // Getting the xPos by using the modulu operator
		index_t yPos = (index_t)(_value/ _geom->Size()[0]) + 1; // Getting the yPos by using division and floor operator
		return multi_index_t(xPos, yPos);
	}
	else return multi_index_t(-1);
}


/// Sets the iterator to the first element
void Iterator::First() {
	_value = 0;
}

/// Goes to the next element of the iterator, disables it if position is end
void Iterator::Next() {
	// Incrementing _value and checking whether it exceeds the maximum of _geom->Size()[0]* _geom->Size()[1]
	_valid = (++_value < _geom->Size()[0]* _geom->Size()[1]) ? true : false;
}


/// Checks if the iterator still has a valid value
bool Iterator::Valid() const {
	return _valid;
}


/// Returns an Iterator that is located left from this one.
// if we are at the left boundary, the cell sees itself
Iterator Iterator::Left() const {
	return (_value % _geom->Size()[0] == 0) ? *this : Iterator(_geom, _value - 1);
}

/// Returns an Iterator that is located right from this one
// If we are at the right boundary, the cell sees itself
Iterator Iterator::Right() const {
	return ((_value+1) % _geom->Size()[0] == 0) ? *this : Iterator(_geom, _value + 1);
}

/// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
Iterator Iterator::Top() const {
	return (_value < _geom->Size()[0]) ? *this : Iterator(_geom, _value - _geom->Size()[0]);
}

/// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
Iterator Iterator::Down() const {
	return (_value >= _geom->Size()[0]*(_geom->Size()[1]-1)) ? *this : Iterator(_geom, _value + _geom->Size()[0]);
}


//------------------------------------------------------------------------------
/** Iterator for interior cells
*/
/// Construct a new InteriorIterator
InteriorIterator::InteriorIterator(const Geometry * geom) {
	_geom = geom;
	First();
	_valid = true;
}


/// Sets the iterator to the first element
void InteriorIterator::First() {
	_value = _geom->Size()[0] + 1;
}

/// Goes to the next element of the iterator, disables it if position is end
void InteriorIterator::Next() {
	// Iterating over inner cells => have to skip a cells when on right border
	// Inner cells:
	// 0 0 0 0 0 0
	// 0 * * * * 0
	// 0 * * * * 0
	// 0 * * * * 0
	// 0 * * * * 0
	// 0 0 0 0 0 0

	_value = ((_value + 2) % _geom->Size()[0] == 0) ? _value + 3 : _value + 1;
	_valid = _value > _geom->Size()[0] * (_geom->Size()[1] - 1) - 2;
}


//------------------------------------------------------------------------------
/** Iterator for domain boundary cells.
*/
/// Constructs a new BoundaryIterator
BoundaryIterator::BoundaryIterator(const Geometry * geom) {
	_geom = geom;
	First();
	_valid = true;
	_boundary = 0;
}


/// Sets the boundary to iterate
void BoundaryIterator::SetBoundary(const index_t & boundary) {
	_boundary = boundary;
}


/// Sets the iterator to the first element
void BoundaryIterator::First() {
	
}

/// Goes to the next element of the iterator, disables it if position is end
void BoundaryIterator::Next() {

}
