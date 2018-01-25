class MGInteriorIterator extends MGIterator {
  
  public MGInteriorIterator(Geometry geom, MultiGrid grid, int cellSize) {
    super(geom, grid, cellSize);
  }
  
  public void next() {
    do {
      _value++;
    } while (_value < _geom.getSize()[0] * _geom.getSize()[1] && (_geom.getFlag()[_value] != ' ' || ((_searchSize == -1 && _grid.getCellSize(_value) != 0) || _grid.getCellSize(_value) != _searchSize)));
    if (_value < _geom.getSize()[0] * _geom.getSize()[1]) _valid = true;
    else _valid = false;
  }  
}