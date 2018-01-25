import java.util.Arrays;

class Grid {
  protected float[] _data;
  protected float[] _offset = new float[2];
  protected Geometry _geom;
  protected float _hInv;
  
  public Grid() {}
  
  public Grid(Geometry geom) {
    _geom = geom;
    _offset[0] = 0.0;
    _offset[1] = 0.0;
    _data = new float[_geom.getSize()[0] * _geom.getSize()[1]];
    _hInv = 1.0 / _geom.getMesh()[0];
  }
  
  public Grid(Geometry geom, float[] offset) {
    _geom = geom;
    _offset = offset;
    _data = new float[_geom.getSize()[0] * _geom.getSize()[1]];
    _hInv = 1.0 / _geom.getMesh()[0];
  }
  
  private Grid(Geometry geom, float[] offset, float[] data) {
    _geom = geom;
    _offset = offset;
    _data = data;
    _hInv = 1.0 / _geom.getMesh()[0];
  }
  
  public Grid copy() {
    return new Grid(_geom, Arrays.copyOf(_offset, 2), Arrays.copyOf(_data, _data.length));
  }
  
  public void printGrid(boolean smooth) {
    for (int i = _geom.getSize()[1] - 1; i >= 0; i--) {
      for (int j = 0; j < _geom.getSize()[0]; j++) {
        if(smooth) print(" " + round(1000*_data[i*_geom.getSize()[0] + j])/1000.0 + " ");
        else print(" " + _data[i*_geom.getSize()[0] + j] + " ");
      }
      println();
    }
    println();
  }
  
  public void initialize(float value) {
    for(int i=0; i<_data.length; i++) _data[i] = value;
  }
  
  float interpolate(float[] pos)  {
    float[] newpos = new float[]{max(0.0, min(pos[0],_geom.getLength()[0])), max(0.0, min(pos[1],_geom.getLength()[1]))};
    
    if(_geom.getFlag()[int(newpos[0] * (_geom.getSize()[0] - 2)/ _geom.getLength()[0] + 1) + int(newpos[1] * (_geom.getSize()[1] - 2)/ _geom.getLength()[1] + 1)*_geom.getSize()[0]] == ' ') {
      float pos_x = newpos[0] * (_geom.getSize()[0] - 2) / _geom.getLength()[0] + 1 - _offset[0];
      float pos_y = newpos[1] * (_geom.getSize()[1] - 2) / _geom.getLength()[1] + 1 - _offset[1];
      int index_x = (int)pos_x;
      int index_y = (int)pos_y;
      
      float val_ll = _data[index_x + index_y*_geom.getSize()[0]];
      float val_lr = _data[index_x + 1 + index_y*_geom.getSize()[0]];
      float val_ul = _data[index_x + (index_y + 1)*_geom.getSize()[0]];
      float val_ur = _data[index_x + 1 + (index_y + 1)*_geom.getSize()[0]];
      
      float prop_x = pos_x - (float)index_x; 
      float prop_y = pos_y - (float)index_y;
      
      return (val_ll*(1.0 - prop_x) + prop_x*val_lr)*(1.0-prop_y) + prop_y*( val_ul*(1.0 - prop_x) + prop_x*val_ur );
    }
    else return 0.0;
  }
  
  public float getCell(Iterator it) {
    return _data[it.getValue()];
  }
  
  public float getCell(int x, int y) {
    return _data[x + y*_geom.getSize()[0]];
  }
  
  public void write2Cell(Iterator it, float value) {
    _data[it.getValue()] = value;
  }
  
  public void write2Cell(int index, float value) {
    _data[index] = value;
  }
  
  public void add2Cell(Iterator it, float value) {
    _data[it.getValue()] += value;
  }
  
  public float dx_l(Iterator it) {
    return (getCell(it) - getCell(it.left())) * _hInv;
  }
  
  public float dx_r(Iterator it) {
    return (getCell(it.right()) - getCell(it)) * _hInv;
  }
  
  public float dy_l(Iterator it) {
    return (getCell(it) - getCell(it.down())) * _hInv;
  }
  
  public float dy_r(Iterator it) {
    return (getCell(it.top()) - getCell(it)) * _hInv;
  }
  
  public float dx_c(Iterator it) {
    return (getCell(it.right()) - getCell(it.down())) * 0.5 * _hInv;
  }
  
  public float dy_c(Iterator it) {
    return (getCell(it.top()) - getCell(it.down())) * 0.5 * _hInv;
  }
  
  public float dxx(Iterator it) {
    return (getCell(it.right()) - 2.0 * getCell(it) + getCell(it.left())) * _hInv * _hInv;
  }

  public float dyy(Iterator it) {
    return (getCell(it.top()) - 2.0 * getCell(it) + getCell(it.down())) * _hInv * _hInv;
  }
  
  public float DC_du2_x(Iterator it, float alpha) {
    float val_u = getCell(it);
    float val_u_r = getCell(it.right());
    float val_u_l = getCell(it.left());
    
    float val_u_cr = (val_u_r + val_u) / 2.0;
    float val_u_cl = (val_u + val_u_l) / 2.0;
    
    return (val_u_cr * val_u_cr - val_u_cl * val_u_cl) * _hInv
        + alpha * (abs(val_u_cr) * (val_u - val_u_r) / 2.0 - abs(val_u_cl) * (val_u_l - val_u) / 2.0) * _hInv;
  }
  
  public float DC_dv2_y(Iterator it, float alpha) {
    float val_v = getCell(it);
    float val_v_t = getCell(it.top());
    float val_v_d = getCell(it.down());
    
    float val_v_ct = (val_v_t + val_v) / 2.0;
    float val_v_cd = (val_v + val_v_d) / 2.0;
    
    return (val_v_ct * val_v_ct - val_v_cd * val_v_cd) * _hInv
        + alpha * (abs(val_v_ct) * (val_v - val_v_t) / 2.0 - abs(val_v_cd) * (val_v_d - val_v) / 2.0) * _hInv;
  }
  
  public float DC_duv_x(Iterator it, float alpha, Grid u) {
    float val_v = getCell(it);
    float val_v_r = getCell(it.right());
    float val_v_l = getCell(it.left());
  
    float val_u_ct = (u.getCell(it.top()) + u.getCell(it)) / 2.0;
    float val_u_ctl = (u.getCell(it.left().top()) + u.getCell(it.left())) / 2.0;
  
    return (val_u_ct * (val_v + val_v_r) / 2.0 - val_u_ctl * (val_v_l + val_v) / 2.0) * _hInv
        + alpha * (abs(val_u_ct) * (val_v - val_v_r) / 2.0 - abs(val_u_ctl) * (val_v_l - val_v) / 2.0) * _hInv;
  }
  
  public float DC_duv_y(Iterator it, float alpha, Grid v) {
    float val_u = getCell(it);
    float val_u_t = getCell(it.top());
    float val_u_d = getCell(it.down());
  
    float val_v_cr = (v.getCell(it.right()) + v.getCell(it)) / 2.0;
    float val_v_cdr = (v.getCell(it.right().down()) + v.getCell(it.down())) / 2.0;
  
    return (val_v_cr * (val_u + val_u_t) / 2.0 - val_v_cdr * (val_u_d + val_u) / 2.0) * _hInv
        + alpha * (abs(val_v_cr) * (val_u - val_u_t) / 2.0 - abs(val_v_cdr) * (val_u_d - val_u) / 2.0) * _hInv;
  }
  
  public float getMax() {
    float max = _data[0];
    for(int i=1; i<_data.length; i++) max = max(max,_data[i]);
    return max;
  }
  
  public float getMin() {
    float min = _data[0];
    for(int i=1; i<_data.length; i++) min = min(min,_data[i]);
    return min;
  }
  
  public float getAbsMax() {
    float max = abs(_data[0]);
    for(int i=1; i<_data.length; i++) max = max(max,abs(_data[i]));
    return max;
  }
  
  public float[] getData() {
    return Arrays.copyOf(_data, _data.length);
  }
  
  public float[] getOffset() {
    return Arrays.copyOf(_offset, 2);
  }
  
  public Geometry getGeometry() {
    return _geom;
  }
}