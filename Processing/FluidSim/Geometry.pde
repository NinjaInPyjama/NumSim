import java.util.Arrays;

class Geometry {
  private int[] _size = new int[2];
  private float[] _length = new float[2];
  private float[] _h = new float[2];
  
  private char[] _flag;
  private float[] _value;
  
  private float[] _velocity = new float[2];
  private float _pressure;
  
  private final String SAVE_FILE_LOCATION;
  
  public Geometry(String filename) {
    SAVE_FILE_LOCATION = "data/" + filename + ".json";
    load();
    
    _value = new float[_size[0]*_size[1]];
    initializeValues();
  }
  
  private void load() {
    JSONObject json = loadJSONObject(SAVE_FILE_LOCATION);
    _size[0] = json.getInt("sizeX");
    _size[1] = json.getInt("sizeY");
    
    _length[0] = json.getFloat("lengthX");
    _length[1] = json.getFloat("lengthY");
    
    _h[0] = _length[0]/_size[0];
    _h[1] = _length[1]/_size[1];
    
    _velocity[0] = json.getFloat("velocityX");
    _velocity[1] = json.getFloat("velocityY");
    
    _pressure = json.getFloat("pressure");
    
    _flag = json.getString("domain").toCharArray();    
  }
  
  void initializeValues() {
    for (int i = 0; i < 4; i++) {
      int start_idx = -1;
      int end_idx = -1;
      int firstID = _size[0] * (_size[1] - 1);
      int lastID = _size[0] * _size[1] - 1;
      int stepID = 1;
      switch (i) {
      case 1: 
        firstID = _size[0]  - 1;
        lastID = _size[0] * _size[1] - 1;
        stepID = _size[0];
        break;
      case 2:
        firstID = 0;
        lastID = _size[0] - 1;
        stepID = 1;
        break;
      case 3:
        firstID = 0;
        lastID = _size[0] * (_size[1] - 1);
        stepID = _size[0];
        break;
      default:
        break;
      }
      for (int j = firstID; j <= lastID; j += stepID) {
        if (start_idx == -1 && (_flag[j] == 'H' || _flag[j] == 'V')) {
          start_idx = j;
        }
        else if (start_idx != -1 && (_flag[j] != 'H' && _flag[j] != 'V')) {
          end_idx = j - stepID;        
          for (int k = start_idx; k <= end_idx; k += stepID) {
            _value[k] = -4.0 / ((end_idx - start_idx + 1.0)*(end_idx - start_idx + 1.0))*_velocity[(i + 1) % 2] * (float(k) - float(start_idx) + 1.0 / 2.0)*(float(k) - float(end_idx) - 1.0 / 2.0);
          }
          start_idx = -1;
          end_idx = -1;
        }
      }
    }
  }
    
  public int[] getSize() {
    return Arrays.copyOf(_size, 2);
  }
  
  public float[] getLength() {
    return Arrays.copyOf(_length, 2);
  }
  
  public float[] getMesh() {
    return Arrays.copyOf(_h, 2);
  }
  
  public char[] getFlag() {
    return _flag;
  }
  
  public float[] getVelocity() {
    return Arrays.copyOf(_velocity, 2);
  }
  
  public float getPressure() {
    return _pressure;
  }
  
  /// Updates the velocity field u
  void updateU(Grid u) {
    // see script, p. 17
    BoundaryIterator bit = new BoundaryIterator(this);
  
    for (bit.first(); bit.isValid(); bit.next()) {
      u.write2Cell(bit, 0.0);
      switch (_flag[bit.getValue()]) {
      case '#': // NOSLIP
        if (_flag[bit.left().getValue()] == ' ') { u.write2Cell(bit.left(), 0.0); u.write2Cell(bit, 0.0); }
        if (_flag[bit.top().getValue()] == ' ') u.write2Cell(bit, -u.getCell(bit.top()));
        if (_flag[bit.down().getValue()] == ' ') u.write2Cell(bit, -u.getCell(bit.down()));
        if (_flag[bit.right().getValue()] == ' ') u.write2Cell(bit, 0.0);
        break;
      case '-': // Horizontal SLIP
        if (_flag[bit.top().getValue()] == ' ') u.write2Cell(bit, u.getCell(bit.top()));
        if (_flag[bit.down().getValue()] == ' ') u.write2Cell(bit, u.getCell(bit.down()));
        if (_flag[bit.left().getValue()] == ' ') u.write2Cell(bit.left(), 0.0);
        if (_flag[bit.right().getValue()] == ' ') u.write2Cell(bit, u.getCell(bit.right()));
        break;
      case '|': // Vertical SLIP
        if (_flag[bit.top().getValue()] == ' ') u.write2Cell(bit, -u.getCell(bit.top()));
        if (_flag[bit.down().getValue()] == ' ') u.write2Cell(bit, -u.getCell(bit.down()));
        if (_flag[bit.left().getValue()] == ' ') u.write2Cell(bit.left(), 0.0);
        if (_flag[bit.right().getValue()] == ' ') u.write2Cell(bit, 0.0);
        break;
      case 'O': // OUTFLOW
        if (_flag[bit.top().getValue()] == ' ') u.write2Cell(bit, u.getCell(bit.top()));
        else if (_flag[bit.down().getValue()] == ' ') u.write2Cell(bit, u.getCell(bit.down()));
        else if (_flag[bit.left().getValue()] == ' ') u.write2Cell(bit, u.getCell(bit.left()));
        else if (_flag[bit.right().getValue()] == ' ') u.write2Cell(bit, u.getCell(bit.right()));
        break;
      case 'V': // Vertical INFLOW
        if (_flag[bit.top().getValue()] == ' ') u.write2Cell(bit, 2.0 * _velocity[0] - u.getCell(bit.top()));
        else if (_flag[bit.down().getValue()] == ' ') u.write2Cell(bit, 2.0 * _velocity[0] - u.getCell(bit.down()));
        else if (_flag[bit.left().getValue()] == ' ') {u.write2Cell(bit, _value[bit.getValue()]); u.write2Cell(bit.left(), _value[bit.getValue()]);}
        else if (_flag[bit.right().getValue()] == ' ') u.write2Cell(bit, _value[bit.getValue()]);
        break;
      case 'H': // Horizontal INFLOW
        if (_flag[bit.top().getValue()] == ' ') u.write2Cell(bit, -u.getCell(bit.top()));
        else if (_flag[bit.down().getValue()] == ' ') u.write2Cell(bit, -u.getCell(bit.down()));
        else if (_flag[bit.left().getValue()] == ' ') u.write2Cell(bit.left(), 0.0);
        else if (_flag[bit.right().getValue()] == ' ') u.write2Cell(bit, 0.0);
        break;
      default:
        break;
      }
    }
  }
  
  /// Updates the velocity field v
  void updateV(Grid v) {
    // see script, p. 17  
    BoundaryIterator bit = new BoundaryIterator(this);
  
    for (bit.first(); bit.isValid(); bit.next()) {
  
      v.write2Cell(bit, 0.0);
      switch (_flag[bit.getValue()]) {
      case '#': // NOSLIP
        if (_flag[bit.down().getValue()] == ' ') { v.write2Cell(bit.down(), 0.0); v.write2Cell(bit, 0.0); }
        if (_flag[bit.left().getValue()] == ' ') v.write2Cell(bit, -v.getCell(bit.left()));
        if (_flag[bit.right().getValue()] == ' ') v.write2Cell(bit, -v.getCell(bit.right()));
        if (_flag[bit.top().getValue()] == ' ') v.write2Cell(bit, 0.0);
        break;
      case '-': // Horizontal SLIP
        if (_flag[bit.left().getValue()] == ' ') v.write2Cell(bit, -v.getCell(bit.left()));
        if (_flag[bit.right().getValue()] == ' ') v.write2Cell(bit, -v.getCell(bit.right()));
        if (_flag[bit.down().getValue()] == ' ') v.write2Cell(bit.down(), 0.0);
        if (_flag[bit.top().getValue()] == ' ') v.write2Cell(bit, 0.0);
        break;
      case '|': // Vertical SLIP
        if (_flag[bit.left().getValue()] == ' ') v.write2Cell(bit, v.getCell(bit.left()));
        if (_flag[bit.right().getValue()] == ' ') v.write2Cell(bit, v.getCell(bit.right()));
        if (_flag[bit.down().getValue()] == ' ') v.write2Cell(bit.down(), 0.0);
        if (_flag[bit.top().getValue()] == ' ') v.write2Cell(bit, 0.0);
        break;
      case 'O': // OUTFLOW
        if (_flag[bit.left().getValue()] == ' ') v.write2Cell(bit, v.getCell(bit.left()));
        else if (_flag[bit.right().getValue()] == ' ') v.write2Cell(bit, v.getCell(bit.right()));
        else if (_flag[bit.down().getValue()] == ' ') v.write2Cell(bit, v.getCell(bit.down()));
        else if (_flag[bit.top().getValue()] == ' ') v.write2Cell(bit, v.getCell(bit.top()));
        break;
      case 'V': // Vertical INFLOW
        if (_flag[bit.left().getValue()] == ' ') v.write2Cell(bit, -v.getCell(bit.left()));
        else if (_flag[bit.right().getValue()] == ' ') v.write2Cell(bit, -v.getCell(bit.right()));
        else if (_flag[bit.down().getValue()] == ' ') v.write2Cell(bit.down(), 0.0);
        else if (_flag[bit.top().getValue()] == ' ') v.write2Cell(bit, 0.0);
        break;
      case 'H': // Horizontal INFLOW
        if (_flag[bit.left().getValue()] == ' ') v.write2Cell(bit, 2.0 * _velocity[1] - v.getCell(bit.left()));
        else if (_flag[bit.right().getValue()] == ' ') v.write2Cell(bit, 2.0 * _velocity[1] - v.getCell(bit.right()));
        else if (_flag[bit.down().getValue()] == ' ') {v.write2Cell(bit, _value[bit.getValue()]); v.write2Cell(bit.down(), _value[bit.getValue()]); }
        else if (_flag[bit.top().getValue()] == ' ') v.write2Cell(bit, _value[bit.getValue()]);
        break;
      default:
        break;
      }
    }
  
  }
  
  /// Updates the pressure field p
  void updateP(Grid p) {
    // see script, p. 20 (p_{0,j} = p{1,j}), ...
    BoundaryIterator bit = new BoundaryIterator(this);
  
  
    for (bit.first(); bit.isValid(); bit.next()) {
      p.write2Cell(bit, 0.0);
      switch (_flag[bit.getValue()]) {
      case '#': // NOSLIP
        if (_flag[bit.top().getValue()] == ' ' || _flag[bit.top().getValue()] == '-') p.write2Cell(bit, p.getCell(bit.top()));
        if (_flag[bit.down().getValue()] == ' ' | _flag[bit.down().getValue()] == '-') {
          if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.down()));
          else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.down())));
        }
        if (_flag[bit.left().getValue()] == ' ') {
          if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.left()));
          else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.left())));
        }
        if (_flag[bit.right().getValue()] == ' ') {
          if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.right()));
          else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.right())));
        }
        break;
      case '-': // Horizontal SLIP
        if (_flag[bit.top().getValue()] == '-' || _flag[bit.down().getValue()] == '-') {
          if (_flag[bit.right().getValue()] == ' ') p.write2Cell(bit, 2.0 * _pressure - p.getCell(bit.right()));
          else p.write2Cell(bit, 2.0 * _pressure - p.getCell(bit.left()));
        }
        else {
          if (_flag[bit.top().getValue()] == ' ') p.write2Cell(bit, p.getCell(bit.top()));
          if (_flag[bit.down().getValue()] == ' ') {
            if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.down()));
            else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.down())));
          }
          if (_flag[bit.left().getValue()] == ' ') {
            if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.left()));
            else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.left())));
          }
          if (_flag[bit.right().getValue()] == ' ') {
            if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.right()));
            else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.right())));
          }
        }
        break;
      case '|': // Vertical SLIP
        if (_flag[bit.right().getValue()] == '|' || _flag[bit.left().getValue()] == '|') {
          if (_flag[bit.top().getValue()] == ' ')p.write2Cell(bit, 2.0 * _pressure - p.getCell(bit.top()));
          else p.write2Cell(bit, 2.0 * _pressure - p.getCell(bit.down()));
        }
        else {
          if (_flag[bit.top().getValue()] == ' ') p.write2Cell(bit, p.getCell(bit.top()));
          if (_flag[bit.down().getValue()] == ' ') {
            if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.down()));
            else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.down())));
          }
          if (_flag[bit.left().getValue()] == ' ') {
            if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.left()));
            else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.left())));
          }
          if (_flag[bit.right().getValue()] == ' ') {
            if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.right()));
            else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.right())));
          }
        }
        break;
      case 'O': // OUTFLOW
        break;
      case 'V': // Vertical INFLOW
        if (_flag[bit.top().getValue()] == ' ') p.write2Cell(bit, p.getCell(bit.top()));
        if (_flag[bit.down().getValue()] == ' ') {
          if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.down()));
          else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.down())));
        }
        if (_flag[bit.left().getValue()] == ' ') {
          if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.left()));
          else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.left())));
        }
        if (_flag[bit.right().getValue()] == ' ') {
          if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.right()));
          else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.right())));
        }
        break;
      case 'H': // Horizontal INFLOW
        if (_flag[bit.top().getValue()] == ' ') p.write2Cell(bit, p.getCell(bit.top()));
        if (_flag[bit.down().getValue()] == ' ') {
          if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.down()));
          else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.down())));
        }
        if (_flag[bit.left().getValue()] == ' ') {
          if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.left()));
          else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.left())));
        }
        if (_flag[bit.right().getValue()] == ' ') {
          if (p.getCell(bit) == 0.0) p.write2Cell(bit, p.getCell(bit.right()));
          else p.write2Cell(bit, 0.5*(p.getCell(bit) + p.getCell(bit.right())));
        }
        break;
      default:
        break;
      }
    }
  }
  
  
  /*
  public void updateU(Grid u) {
    BoundaryIterator bit = new BoundaryIterator(this);
    
    bit.setBoundary(1);
    for(bit.first(); bit.isValid(); bit.next()) {
      u.write2Cell(bit.left(), u.getCell(bit.left().left()));
      u.write2Cell(bit,  u.getCell(bit.left()));
    }
    
    bit.setBoundary(3);
    for(bit.first(); bit.isValid(); bit.next()) {
      u.write2Cell(bit,  u.getCell(bit.right()));
    }
    
    bit.setBoundary(0);
    for(bit.first(); bit.isValid(); bit.next()) {
      u.write2Cell(bit, - u.getCell(bit.down()));
    }
  
    bit.setBoundary(2);
    for(bit.first(); bit.isValid(); bit.next()) {
      u.write2Cell(bit, - u.getCell(bit.top()));
    }
  }
  
  public void updateV(Grid v) {
    BoundaryIterator bit = new BoundaryIterator(this);
    
    bit.setBoundary(1);
    for(bit.first(); bit.isValid(); bit.next()) {
      v.write2Cell(bit,  v.getCell(bit.left()));
    }
    
    bit.setBoundary(3);
    for(bit.first(); bit.isValid(); bit.next()) {
      v.write2Cell(bit, 0.0);
    }
    
    bit.setBoundary(0);
    for(bit.first(); bit.isValid(); bit.next()) {
      v.write2Cell(bit, 0.0);
      v.write2Cell(bit.down(), 0.0);
    }
    
    bit.setBoundary(2);
    for(bit.first(); bit.isValid(); bit.next()) {
      v.write2Cell(bit, 0.0);
    }
  }
  
  public void updateP(Grid p) {
    BoundaryIterator bit = new BoundaryIterator(this);
    
    bit.setBoundary(1);
    for(bit.first(); bit.isValid(); bit.next()) {
      p.write2Cell(bit, 0.0);
    }
    
    bit.setBoundary(3);
    for(bit.first(); bit.isValid(); bit.next()) {
      p.write2Cell(bit, 0.1);
    }
    
    bit.setBoundary(0);
    for(bit.first(); bit.isValid(); bit.next()) {
      p.write2Cell(bit, p.getCell(bit.down()));
    }
    
    bit.setBoundary(2);
    for(bit.first(); bit.isValid(); bit.next()) {
      p.write2Cell(bit, p.getCell(bit.top()));
    }
  }
  
  /*
  public void updateU(Grid u) {
    BoundaryIterator bit = new BoundaryIterator(this);
    
    bit.setBoundary(1);
    for(bit.first(); bit.isValid(); bit.next()) {
      u.write2Cell(bit, 0.0);
      u.write2Cell(bit.left(), 0.0);
    }
    
    bit.setBoundary(3);
    for(bit.first(); bit.isValid(); bit.next()) {
      u.write2Cell(bit, 0.0);
    }
    
    bit.setBoundary(0);
    for(bit.first(); bit.isValid(); bit.next()) {
      u.write2Cell(bit, 2.0 - u.getCell(bit.down()));
    }
  
    bit.setBoundary(2);
    for(bit.first(); bit.isValid(); bit.next()) {
      u.write2Cell(bit, - u.getCell(bit.top()));
    }
  }
  
  public void updateV(Grid v) {
    BoundaryIterator bit = new BoundaryIterator(this);
    
    bit.setBoundary(0);
    for(bit.first(); bit.isValid(); bit.next()) {
      v.write2Cell(bit, 0.0);
      v.write2Cell(bit.down(), 0.0);
    }
    
    bit.setBoundary(2);
    for(bit.first(); bit.isValid(); bit.next()) {
      v.write2Cell(bit, 0.0);
    }
    
    bit.setBoundary(3);
    for(bit.first(); bit.isValid(); bit.next()) {
      v.write2Cell(bit, - v.getCell(bit.right()));
    }
    
    bit.setBoundary(1);
    for(bit.first(); bit.isValid(); bit.next()) {
      v.write2Cell(bit, - v.getCell(bit.left()));
    }
  }
  
  public void updateP(Grid p) {
    BoundaryIterator bit = new BoundaryIterator(this);
    
    bit.setBoundary(0);
    for(bit.first(); bit.isValid(); bit.next()) {
      p.write2Cell(bit, p.getCell(bit.down()));
    }
    
    bit.setBoundary(1);
    for(bit.first(); bit.isValid(); bit.next()) {
      p.write2Cell(bit, p.getCell(bit.left()));
    }
    
    bit.setBoundary(2);
    for(bit.first(); bit.isValid(); bit.next()) {
      p.write2Cell(bit, p.getCell(bit.top()));
    }
    
    bit.setBoundary(3);
    for(bit.first(); bit.isValid(); bit.next()) {
      p.write2Cell(bit, p.getCell(bit.right()));
    }
  }
  */
}