class MGIterator extends Iterator {
  protected Geometry _geom;
  protected MultiGrid _grid;
  protected int _searchSize;
  
  public MGIterator(){}
  
  public MGIterator(Geometry geom, MultiGrid grid, int cellSize) {
    _geom = geom;
    _grid = grid;
    _searchSize = cellSize;
    first();
  }
  
  public MGIterator(Geometry geom, MultiGrid grid, int cellSize, int value) {
    _geom = geom;
    _grid = grid;
    _searchSize = cellSize;
    _value = value;
    _valid = true;
  }
  
  public void first() {
    _value = -1;
    _valid = true;
    next();
  }
  
  public MGIterator left() {
    int currCellSize = _grid.getCellSize(_value);
    if(_value % _geom.getSize()[0] == 0) return this;
    else {
      if(_grid.getCellSize(_value-1) != 0) return new MGIterator(_geom, _grid, -1, _value-1);
      else if(currCellSize != 1 && _grid.getCellSize(_value-int(0.5*currCellSize)) != 0) return new MGIterator(_geom, _grid, -1, _value-int(0.5*currCellSize));
      else if(_grid.getCellSize(_value-currCellSize) != 0) return new MGIterator(_geom, _grid, -1, _value-currCellSize);
      else if(_grid.getCellSize(_value-2*currCellSize) != 0) return new MGIterator(_geom, _grid, -1, _value-2*currCellSize);
      else if(_grid.getCellSize(_value-currCellSize*(2+_geom.getSize()[0])) != 0) return new MGIterator(_geom, _grid, -1, _value-currCellSize*(2+_geom.getSize()[0]));
      else return this;
    }
  }
    
  public MGIterator right() {
    int currCellSize = _grid.getCellSize(_value);
    if(_value+1 % _geom.getSize()[0] == 0) return this;
    else {
      if(_grid.getCellSize(_value+1) != 0) return new MGIterator(_geom, _grid, -1, _value+1);
      else if(_grid.getCellSize(_value+currCellSize) != 0) return new MGIterator(_geom, _grid, -1, _value+currCellSize);
      else if(_grid.getCellSize(_value+currCellSize*(1-_geom.getSize()[0])) != 0) return new MGIterator(_geom, _grid, -1, _value+currCellSize*(1-_geom.getSize()[0]));
      else return this;
    }
  }
  
  public MGIterator top() {
    int currCellSize = _grid.getCellSize(_value);
    if(_value >= _geom.getSize()[0] * (_geom.getSize()[1] - 1)) return this;
    else {
      if(_grid.getCellSize(_value+_geom.getSize()[0]) != 0) return new MGIterator(_geom, _grid, -1, _value+_geom.getSize()[0]);
      else if(_grid.getCellSize(_value+currCellSize*_geom.getSize()[0]) != 0) return new MGIterator(_geom, _grid, -1, _value+currCellSize*_geom.getSize()[0]);
      else if(_grid.getCellSize(_value+currCellSize*(_geom.getSize()[0]-1)) != 0) return new MGIterator(_geom, _grid, -1, _value+currCellSize*(_geom.getSize()[0]-1));
      else return this;
    }
  }
  
  public MGIterator down() {
    int currCellSize = _grid.getCellSize(_value);
    if(_value < _geom.getSize()[0]) return this;
    else {
      if(_grid.getCellSize(_value-_geom.getSize()[0]) != 0) return new MGIterator(_geom, _grid, -1, _value-_geom.getSize()[0]);
      else if(currCellSize != 1 && _grid.getCellSize(_value-int(0.5*currCellSize*_geom.getSize()[0])) != 0) return new MGIterator(_geom, _grid, -1, _value-int(0.5*currCellSize*_geom.getSize()[0]));
      else if(_grid.getCellSize(_value-currCellSize*_geom.getSize()[0]) != 0) return new MGIterator(_geom, _grid, -1, _value-currCellSize*_geom.getSize()[0]);
      else if(_grid.getCellSize(_value-2*currCellSize*_geom.getSize()[0]) != 0) return new MGIterator(_geom, _grid, -1, _value-2*currCellSize*_geom.getSize()[0]);
      else if(_grid.getCellSize(_value-currCellSize*(2*_geom.getSize()[0]+1)) != 0) return new MGIterator(_geom, _grid, -1, _value-currCellSize*(2*_geom.getSize()[0]+1));
      else return this;
    }
  }
}