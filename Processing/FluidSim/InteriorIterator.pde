class InteriorIterator extends Iterator {
  
  public InteriorIterator(Geometry geom){
    _geom = geom;
    _valid = true;
    first();
  }
  
  public void first() {
    _value = -1;
    next();
  }
  
  public void next() {
    do {
      _value++;
    } while (_value < _geom.getSize()[0] * _geom.getSize()[1] && _geom.getFlag()[_value] != ' ');
    if (_value < _geom.getSize()[0] * _geom.getSize()[1]) _valid = true;
    else _valid = false;
  }
}