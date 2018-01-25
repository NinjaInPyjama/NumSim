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
/** Iterator for multigrid algorithm
*/
/// Construct a new InteriorIterator
MGIterator::MGIterator(const Geometry *geom, const MultiGrid *multigrid, const index_t cellsize){
    _geom = geom;
    _valid = true;
    _multigrid = multigrid;
    _cellsize = cellsize;
    First();   
}

/// Constructs a new BoundaryIterator to start iterating from chosen value
MGIterator::MGIterator(const Geometry *geom, const MultiGrid *multigrid, const index_t cellsize, const index_t &value){
    _geom = geom;
    _valid = true;
    _multigrid = multigrid;
    _cellsize = cellsize;
    _value = value;
}

/// Sets the iterator to the first element
void MGInteriorIterator::First(){
 	_value = -1;
	Next();   
}

/// Returns an Iterator that is located left from this one.
// if we are at the left boundary, the cell sees itself
// if there is more than one cell left, we choose the lower one
MGIterator MGIterator::Left() const{
    index_t current_size = _multigrid->CellSize(value);
    if (_value % _geom->Size()[0] == 0) return *this;
    else if (_multigrid->CellSize(value-1) != 0) return MGIterator(_geom, _multigrid, _multigrid-CellSize(value-1), value-1);
    else if (current_size!=1 && _multigrid->CellSize(value-current_size/2) != 0) return MGIterator(_geom, _multigrid, _multigrid->CellSize(value-current_size/2), value-current_size/2);
    else if (_multigrid->CellSize(value-current_size) != 0) return MGIterator(_geom, _multigrid, _multigrid->CellSize(value-current_size), value-current_size);
    else if (_multigrid->CellSize(value-2*current_size) != 0) return MGIterator(_geom, _multigrid, _multigrid->CellSize(value-2*current_size), value-2*current_size);
    else if (_multigrid->CellSize(value-2*current_size-current_size*_geom->Size()[0]) != 0) return MGIterator(_geom, _multigrid, _multigrid->CellSize(value-2*current_size-current_size*_geom->Size()[0]), value-2*current_size-current_size*_geom->Size()[0]);
    else return *this;
}

/// Returns an Iterator that is located right from this one
// If we are at the right boundary, the cell sees itself
// if there is more than one cell right, we choose the lower one
MGIterator MGIterator::Right() const{
    index_t current_size = _multigrid->CellSize(value);
    if ((_value + 1)% _geom->Size()[0] == 0) return *this;
    else if (_multigrid->CellSize(value + current_size) != 0) return MGIterator(_geom, _multigrid, _multigrid-CellSize(value+current_size), value+current_size);
    else if (_multigrid->CellSize(value + current_size - current_size*_geom->Size()[0]) != 0) return MGIterator(_geom, _multigrid, _multigrid->CellSize(value + current_size - current_size*_geom->Size()[0]), value + current_size - current_size*_geom->Size()[0]);
    else return *this;
}

/// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
// if there is more than one cell at the top, we choose the left one
MGIterator MGIterator::Top() const{
    index_t current_size = _multigrid->CellSize(value);
    if (_value >= _geom->Size()[0]*(_geom->Size()[1]-1)) return *this;
    else if (_multigrid->CellSize(value + current_size*_geom->Size()[0]) != 0) return MGIterator(_geom, _multigrid, _multigrid->CellSize(value + current_size*_geom->Size()[0]), value + current_size*_geom->Size()[0]);
    else if (_multigrid->CellSize(value + current_size*_geom->Size()[0] - current_size) != 0) return MGIterator(_geom, _multigrid, _multigrid->CellSize(value + current_size*_geom->Size()[0] - current_size), value + current_size*_geom->Size()[0] - current_size);
    else return *this;
}

/// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
// if there is more than one cell below, we choose the left one
MGIterator MGIterator::Down() const{
    index_t current_size = _multigrid->CellSize(value);
    if (_value < _geom->Size()[0]) return *this;
    else if (_multigrid->CellSize(value-_geom->Size()[0]) != 0) return MGIterator(_geom, _multigrid, _multigrid-CellSize(value-_geom->Size()[0])), value-_geom->Size()[0]));
    else if (current_size!=1 && _multigrid->CellSize(value-_geom->Size()[0]*current_size/2) != 0) return MGIterator(_geom, _multigrid, _multigrid->CellSize(value-_geom->Size()[0]*current_size/2), value-_geom->Size()[0]*current_size/2);
    else if (_multigrid->CellSize(value-_geom->Size()[0]*current_size) != 0) return MGIterator(_geom, _multigrid, _multigrid->CellSize(value-_geom->Size()[0]*current_size), value-_geom->Size()[0]*current_size);
    else if (_multigrid->CellSize(value-2*_geom->Size()[0]*current_size) != 0) return MGIterator(_geom, _multigrid, _multigrid->CellSize(value-2*_geom->Size()[0]*current_size), value-2*_geom->Size()[0]*current_size);
    else if (_multigrid->CellSize(value-current_size-2*current_size*_geom->Size()[0]) != 0) return MGIterator(_geom, _multigrid, _multigrid->CellSize(value-current_size-2*current_size*_geom->Size()[0]), value-current_size-2*current_size*_geom->Size()[0]);
    else return *this;
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
/// Goes to the next element of the iterator, disables it if position is end
void MGInteriorIterator::Next(){
 	do {
		_value++;
	} while (_value < _geom->Size()[0] * _geom->Size()[1] && (_geom->Flag()[_value] != ' ' || ((_cellsize == -1 && _multigrid->CellSize(this) != 0) || _multigrid->CellSize(this) != _cellsize)));
	if (_value < _geom->Size()[0] * _geom->Size()[1]) _valid = true;
	else _valid = false;   
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
    _boundary = 0;
    First();
}

/// Sets boundary
void MGBoundaryIterator::Boundary(index_t boundary){
    _boundary = boundary;
}

/// Sets the iterator to the first element
void MGBoundaryIterator::First(){
    switch(_boundary) {
      case 0: // Right Boundary
        _value = 2*_geom->Size()[0] - 1;
        break;
      case 1: // Bottom Boundary
        _value = 1;
        break;
      case 2: // Left Boundary
      _value = _geom->Size()[0];
        break;
      case 3: // Top Boundary
        _value = _geom->Size()[0]*(_geom->Size()[1] - 1) + 1;
        break; 
      case 4: // Inner Boundaries
        _value = _geom->Size()[0] + 1;
        break;
      default:
        break;
    }
    _valid = true;
}

/// Goes to the next element of the iterator, disables it if position is end
void MGBoundaryIterator::Next(){
    	switch(_boundary) {
      case 0: // Right Boundary
        do {
          _value += _geom->Size()[0];
        } while (_value < _geom->Size()[0] * (_geom->Size()[1]-1) && ((_cellsize == -1 && _grid->CellSize(_value) != 0) || _grid->CellSize(_value) != _cellsize));
        if (_value < _geom->Size()[0] * (_geom->Size()[1]-1)) _valid = true;
        else _valid = false;
        break;
      case 1: // Bottom Boundary
        do {
          _value++;
        } while (_value < _geom->Size()[0] - 1 && ((_cellsize == -1 && _grid->CellSize(_value) != 0) || _grid->CellSize(_value) != _cellsize));
        if (_value < _geom->Size()[0] - 1) _valid = true;
        else _valid = false;
        break;
      case 2: // Left Boundary
        do {
          _value += _geom->Size()[0];
        } while (_value < _geom->Size()[0] * (_geom->Size()[1]-1) && (_geom->Flag()[_value] == ' ' || ((_cellsize == -1 && _grid->CellSize(_value) != 0) || _grid->CellSize(_value) != _cellsize)));
        if (_value < _geom->Size()[0] * (_geom->Size()[1]-1)) _valid = true;
        else _valid = false;
        break;
      case 3: // Top Boundary
        do {
          _value++;
        } while (_value < _geom->Size()[0] * _geom->Size()[1] - 1 && (_geom->Flag()[_value] == ' ' || ((_cellsize == -1 && _grid->CellSize(_value) != 0) || _grid->CellSize(_value) != _cellsize)));
        if (_value < _geom->Size()[0] * _geom->Size()[1] - 1) _valid = true;
        else _valid = false;
        break; 
      case 4: // Inner Boundaries
        do {
          _value++;
        } while (_value < _geom->Size()[0] * (_geom->Size()[1]-1) - 1 && (_value%_geom->Size()[0] == 0 || (_value-1)%_geom->Size()[0] == 0 || _geom->Flag()[_value] == ' ' || ((_cellsize == -1 && _grid->CellSize(_value) != 0) || _grid->CellSize(_value) != _cellsize)));
        if (_value < _geom->Size()[0] * (_geom->Size()[1]-1) - 1) _valid = true;
        else _valid = false;
        break;
      default:
        break;
    }
}