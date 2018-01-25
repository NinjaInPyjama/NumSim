class Iterator {
  protected Geometry _geom;
  protected int _value;
  protected boolean _valid;
  
  public Iterator(){}
  
  public Iterator(Geometry geom){
    _geom = geom;
    first();
  }
  
   public Iterator(Geometry geom, int value){
    _geom = geom;
    _value = value;
    _valid = true;
  }
  
  public int[] pos() {
    if(_valid) {
      int xPos = (_value % _geom.getSize()[0]);
      int yPos = (int)(float(_value) / float(_geom.getSize()[0]));
      return new int[]{xPos, yPos};
    }
    else return new int[]{-1, -1};
  }
  
  public void first() {
    _value = 0;
    _valid = true;
  }
  
  public void next() {
    _value++;
    _valid = _value < _geom.getSize()[0] * _geom.getSize()[1];
  }
  
  public int getValue() {
    return _value;
  }
  
  public boolean isValid() {
    return _valid;
  }
  
  public Iterator left() {
    if(_value % _geom.getSize()[0] == 0) return this;
    else return new Iterator(_geom,_value - 1);
  }
  
  public Iterator right() {
    if((_value + 1) % _geom.getSize()[0] == 0) return this;
    else return new Iterator(_geom,_value + 1);
  }
  
  public Iterator top() {
    if(_value >= _geom.getSize()[0] * (_geom.getSize()[1] - 1)) return this;
    else return new Iterator(_geom,_value + _geom.getSize()[0]);
  }
  
  public Iterator down() {
    if(_value < _geom.getSize()[0]) return this;
    else return new Iterator(_geom,_value - _geom.getSize()[0]);
  }
  
}