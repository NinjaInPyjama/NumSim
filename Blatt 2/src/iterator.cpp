#include "iterator.hpp"

/// Default constructor
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
		index_t xPos = (_value % _geom->Size()[0]) ; // Getting the xPos by using the modulu operator
		index_t yPos = (index_t)(real_t(_value) / real_t(_geom->Size()[0])) ; // Getting the yPos by using division and floor operator
		return multi_index_t(xPos, yPos);
	}
	else return multi_index_t(-1);
}


/// Sets the iterator to the first element
void Iterator::First() {
	_value = 0;
	_valid = true;
}

/// Goes to the next element of the iterator, disables it if position is end
void Iterator::Next() {
	// Incrementing _value and checking whether it exceeds the maximum of _geom->Size()[0] * _geom->Size()[1]
	_valid = (++_value < _geom->Size()[0] * _geom->Size()[1]) ? true : false;
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
	return ((_value + 1) % _geom->Size()[0] == 0) ? *this : Iterator(_geom, _value + 1);
}

/// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
Iterator Iterator::Top() const {
	return (_value >= _geom->Size()[0] * (_geom->Size()[1] - 1)) ? *this : Iterator(_geom, _value + _geom->Size()[0]);
}

/// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
Iterator Iterator::Down() const {
	return (_value < _geom->Size()[0]) ? *this : Iterator(_geom, _value - _geom->Size()[0]);
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
	_valid = true;
}

/// Goes to the next element of the iterator, disables it if position is end
// Iterating over inner cells => has to skip cells when on right border
// Inner cells (*):
// 0 0 0 0 0 0
// 0 * * * * 0
// 0 * * * * 0
// 0 * * * * 0
// 0 * * * * 0
// 0 0 0 0 0 0
void InteriorIterator::Next() {
	_value = ((_value + 2) % _geom->Size()[0] == 0) ? _value + 3 : _value + 1;
	_valid = _value <= _geom->Size()[0] * (_geom->Size()[1] - 1) - 2;
}




//------------------------------------------------------------------------------
/** Iterator for interior cells for the red black solver
*/

/// Construct a new RedBlackIterator
RedBlackIterator::RedBlackIterator(const Geometry * geom, bool rb) {
	_rb = rb; //true means red cycle
	_geom = geom;
	First();
	_valid = true;
}


/// Sets the iterator to the first red element
void RedBlackIterator::First() {
	_value = (_geom->RedBlack() == _rb) ? _geom->Size()[0] + 1 : _geom->Size()[0] + 2;
	_valid = true;
}

/// Goes to the next element of the iterator, disables it if position is end
// Iterating over only red or black cells
void RedBlackIterator::Next() {
    if (_geom->Size()[0] % 2 == 0) {
		if ((_value + 3) % _geom->Size()[0] == 0) {
			_value = _value + 3;
		}
		else if ((_value + 2) % _geom->Size()[0] == 0) {
			_value = _value + 1;
		}
	}
	else {
		if ((_value + 3) % _geom->Size()[0] == 0 || (_value + 2) % _geom->Size()[0] == 0) {
			_value = _value + 2;
		}
	}
	_value = _value + 2;
    _valid = _value < _geom->Size()[0] * (_geom->Size()[1] - 1) - 1;
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
// boundary = 0 : iterating over upper boundary
// boundary = 1 : iterating over right boundary
// boundary = 2 : iterating over lower boundary
// boundary = 3 : iterating over left boundary
// every other input will later be treated as 0.
void BoundaryIterator::SetBoundary(const index_t & boundary) {
	_boundary = boundary;
}


/// Sets the iterator to the first element
void BoundaryIterator::First() {
	_valid = true;
	switch(_boundary) {
        case 0: // upper boundary
			_value = _geom->Size()[0]*(_geom->Size()[1] - 1);
			break;
		case 1: // right boundary
			_value = _geom->Size()[0] - 1;
			break;
		case 2: // lower boundary
			_value = 0;
			break;
		case 3: // left boundary
			_value = 0;
			break;
		default: // simulates top boundary
			_value = _geom->Size()[0]*(_geom->Size()[1] - 1);
			break;
	}
}

/// Goes to the next element of the iterator, disables it if position is end
// Iterating over boundary cells without corners
// Boundary cells (*):
// 0 * * * * 0
// * 0 0 0 0 *
// * 0 0 0 0 *
// * 0 0 0 0 *
// * 0 0 0 0 *
// 0 * * * * 0
void BoundaryIterator::Next() {
	switch(_boundary) {
		case 0: // upper boundary
			_value++;
			_valid = (_value < _geom->Size()[0] * _geom->Size()[1]) ? true : false;		
			break;
		case 1: // right boundary
			_value = _value + _geom->Size()[0];
			_valid = (_value < _geom->Size()[0] * _geom->Size()[1]) ? true : false;
			break;
		case 2: // lower boundary
			_value++;
			_valid = (_value < _geom->Size()[0]) ? true : false;
			break;
		case 3: // left boundary
			_value = _value + _geom->Size()[0];
			_valid = (_value < _geom->Size()[0] * (_geom->Size()[1] - 1) + 1) ? true : false;		
			break;
		default: // simulates upper boundary
			_value++;
			_valid = (_value < _geom->Size()[0] * _geom->Size()[1]) ? true : false;
			break;
	}
}
