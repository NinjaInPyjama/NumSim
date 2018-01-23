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
	return ((_value+1) % _geom->Size()[0] == 0) ? *this : Iterator(_geom, _value + 1);
}

/// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
Iterator Iterator::Top() const {
	return (_value >= _geom->Size()[0]*(_geom->Size()[1]-1)) ? *this : Iterator(_geom, _value + _geom->Size()[0]);
}

/// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
Iterator Iterator::Down() const {
	return (_value < _geom->Size()[0]) ? *this : Iterator(_geom, _value - _geom->Size()[0]);
}


//------------------------------------------------------------------------------
/** Iterator for all cells running backwards linewise
*/

/// Construct a new BackwardsIterator
BackwardsIterator::BackwardsIterator(const Geometry * geom) {
	_geom = geom;
	_valid = true;
	First();
}


/// Sets the iterator to the last element of the first line
void BackwardsIterator::First() {
	_value = _geom->Size()[0]-1;
}

/// Goes to the next element of the iterator, disables it if position is end
// Iterating over inner cells
void BackwardsIterator::Next() {
	
    if(_value % _geom->Size()[0] == 0) {
        _value += 2*_geom->Size()[0] - 1;
    }
    else {
        _value--;
    }
    if(_value > _geom->Size()[0]*_geom->Size()[1]-1) _valid = false;
}



//------------------------------------------------------------------------------
/** Iterator for interior cells
*/

/// Construct a new InteriorIterator
InteriorIterator::InteriorIterator(const Geometry * geom) {
	_geom = geom;
	_valid = true;
	First();
}


/// Sets the iterator to the first element
void InteriorIterator::First() {
	_value = -1;
	Next();
}

/// Goes to the next element of the iterator, disables it if position is end
// Iterating over inner cells
void InteriorIterator::Next() {
	do {
		_value++;
	} while (_geom->Flag()[_value] != ' ' && _value < _geom->Size()[0] * _geom->Size()[1]);
	if (_value < _geom->Size()[0] * _geom->Size()[1]) _valid = true;
	else _valid = false;
}

//------------------------------------------------------------------------------
/** Iterator for domain boundary cells.
*/
/// Constructs a new BoundaryIterator
BoundaryIterator::BoundaryIterator(const Geometry * geom) {
	_geom = geom;
	_valid = true;
	First();
}


/// Sets the iterator to the first element
void BoundaryIterator::First() {
	_value = -1;
	Next();
}

/// Goes to the next element of the iterator, disables it if position is end
// Iterating over boundary cells 
void BoundaryIterator::Next() {
	do {
		_value++;
	} while (_geom->Flag()[_value] == ' ' && _value < _geom->Size()[0] * _geom->Size()[1]);
	if (_value < _geom->Size()[0] * _geom->Size()[1]) _valid = true;
	else _valid = false;
}

//------------------------------------------------------------------------------
/** Iterator for interior cells in multigrid algorithm
*/
/// Construct a new InteriorIterator
MGInteriorIterator::MGInteriorIterator(const Geometry *geom, const MultiGrid *multigrid, index_t cellsize){
 	_geom = geom;
	_valid = true;
        _multigrid = multigrid;
        _cellsize = cellsize;
	First();   
}

/// Constructs a new BoundaryIterator to start iterating from chosen value
MGInteriorIterator::MGInteriorIterator(const Geometry *geom, const MultiGrid *multigrid, index_t cellsize, index_t value){
    
}

/// Sets the iterator to the first element
void MGInteriorIterator::First(){
 	_value = -1;
	Next();   
}
/// Goes to the next element of the iterator, disables it if position is end
void MGInteriorIterator::Next(){
 	do {
		_value++;
	} while (((_cellsize == -1 && _multigrid->CellSize(this) != 0) || _multigrid->CellSize(this) != _cellsize) && _geom->Flag()[_value] != ' ' && _value < _geom->Size()[0] * _geom->Size()[1]);
	if (_value < _geom->Size()[0] * _geom->Size()[1]) _valid = true;
	else _valid = false;   
}
/// Returns an Iterator that is located left from this one.
// if we are at the left boundary, the cell sees itself
// if there is more than one cell left, we choose the lower one
MGInteriorIterator MGInteriorIterator::Left() const{
    
}

/// Returns an Iterator that is located right from this one
// If we are at the right boundary, the cell sees itself
// if there is more than one cell right, we choose the lower one
MGInteriorIterator MGInteriorIterator::Right() const{
    
}

/// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
// if there is more than one cell at the top, we choose the left one
MGInteriorIterator MGInteriorIterator::Top() const{
    
}

/// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
// if there is more than one cell below, we choose the left one
MGInteriorIterator MGInteriorIterator::Down() const{
    
}

//------------------------------------------------------------------------------
/** Iterator for domain boundary cells in multigrid algorithm
*/
/// Constructs a new BoundaryIterator
MGBoundaryIterator::MGBoundaryIterator(const Geometry *geom, const MultiGrid *multigrid, index_t cellsize){
    	_geom = geom;
	_valid = true;
        _multigrid = multigrid;
        _cellsize = cellsize;
	First();
}

/// Constructs a new BoundaryIterator to start iterating from chosen value
MGBoundaryIterator::MGBoundaryIterator(const Geometry *geom, const MultiGrid *multigrid, index_t cellsize, index_t value){
    
}

/// Sets the iterator to the first element
void MGBoundaryIterator::First(){
    	_value = -1;
	Next();
}

/// Goes to the next element of the iterator, disables it if position is end
void MGBoundaryIterator::Next(){
    	do {
		_value++;
	} while (((_cellsize == -1 && _multigrid->CellSize(this) != 0) || _multigrid->CellSize(this) != _cellsize) && _geom->Flag()[_value] == ' ' && _value < _geom->Size()[0] * _geom->Size()[1]);
	if (_value < _geom->Size()[0] * _geom->Size()[1]) _valid = true;
	else _valid = false;
}

/// Returns an Iterator that is located left from this one.
// if we are at the left boundary, the cell sees itself
// if there is more than one cell left, we choose the lower one
MGBoundaryIterator MGBoundaryIterator::Left() const{
    
}

/// Returns an Iterator that is located right from this one
// If we are at the right boundary, the cell sees itself
// if there is more than one cell right, we choose the lower one
MGBoundaryIterator MGBoundaryIterator::Right() const{
    
}

/// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
// if there is more than one cell at the top, we choose the left one
MGBoundaryIterator MGBoundaryIterator::Top() const{
    
}

/// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
// if there is more than one cell below, we choose the left one
MGBoundaryIterator MGBoundaryIterator::Down() const{
    
}
