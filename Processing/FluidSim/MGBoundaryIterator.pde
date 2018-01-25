class MGBoundaryIterator extends MGIterator {
  
  private int _boundary;

  public MGBoundaryIterator(Geometry geom, MultiGrid grid, int cellSize) {
    super(geom, grid, cellSize);
    _boundary = 0;
  }
  
  public void setBoundary(int boundary) {
    _boundary = boundary;
  }
  
  public void first() {
    switch(_boundary) {
      case 0: // Right Boundary
        _value = 2*_geom.getSize()[0] - 1;
        break;
      case 1: // Bottom Boundary
        _value = 1;
        break;
      case 2: // Left Boundary
      _value = _geom.getSize()[0];
        break;
      case 3: // Top Boundary
        _value = _geom.getSize()[0]*(_geom.getSize()[1] - 1) + 1;
        break; 
      case 4: // Inner Boundaries
        _value = _geom.getSize()[0] + 1;
        break;
      default:
        break;
    }
    _valid = true;
  }
  
  public void next() {
    switch(_boundary) {
      case 0: // Right Boundary
        do {
          _value += _geom.getSize()[0];
        } while (_value < _geom.getSize()[0] * (_geom.getSize()[1]-1) && ((_searchSize == -1 && _grid.getCellSize(_value) != 0) || _grid.getCellSize(_value) != _searchSize));
        if (_value < _geom.getSize()[0] * (_geom.getSize()[1]-1)) _valid = true;
        else _valid = false;
        break;
      case 1: // Bottom Boundary
        do {
          _value++;
        } while (_value < _geom.getSize()[0] - 1 && ((_searchSize == -1 && _grid.getCellSize(_value) != 0) || _grid.getCellSize(_value) != _searchSize));
        if (_value < _geom.getSize()[0] - 1) _valid = true;
        else _valid = false;
        break;
      case 2: // Left Boundary
        do {
          _value += _geom.getSize()[0];
        } while (_value < _geom.getSize()[0] * (_geom.getSize()[1]-1) && (_geom.getFlag()[_value] == ' ' || ((_searchSize == -1 && _grid.getCellSize(_value) != 0) || _grid.getCellSize(_value) != _searchSize)));
        if (_value < _geom.getSize()[0] * (_geom.getSize()[1]-1)) _valid = true;
        else _valid = false;
        break;
      case 3: // Top Boundary
        do {
          _value++;
        } while (_value < _geom.getSize()[0] * _geom.getSize()[1] - 1 && (_geom.getFlag()[_value] == ' ' || ((_searchSize == -1 && _grid.getCellSize(_value) != 0) || _grid.getCellSize(_value) != _searchSize)));
        if (_value < _geom.getSize()[0] * _geom.getSize()[1] - 1) _valid = true;
        else _valid = false;
        break; 
      case 4: // Inner Boundaries
        do {
          _value++;
        } while (_value < _geom.getSize()[0] * (_geom.getSize()[1]-1) - 1 && (_value%_geom.getSize()[0] == 0 || (_value-1)%_geom.getSize()[0] == 0 || _geom.getFlag()[_value] == ' ' || ((_searchSize == -1 && _grid.getCellSize(_value) != 0) || _grid.getCellSize(_value) != _searchSize)));
        if (_value < _geom.getSize()[0] * (_geom.getSize()[1]-1) - 1) _valid = true;
        else _valid = false;
        break;
      default:
        break;
    }
  }
}