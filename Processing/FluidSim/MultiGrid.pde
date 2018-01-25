class MultiGrid extends Grid{
  private int[] _cellSize;
  
  public MultiGrid(Geometry geom) {
    _geom = geom;
    _cellSize = new int[_geom.getSize()[0] * _geom.getSize()[1]];
    for(int i=0; i<_cellSize.length; i++) _cellSize[i] = 1;
    _data = new float[_geom.getSize()[0] * _geom.getSize()[1]];
    initialize(0.0);
  }
  
  public MultiGrid(Geometry geom, float[] data) {
    _geom = geom;
    _data = data;
    _cellSize = new int[_geom.getSize()[0] * _geom.getSize()[1]];
    for(int i=0; i<_cellSize.length; i++) _cellSize[i] = 1;
  }
  
  public MultiGrid(Grid grid) {
    _geom = grid._geom;
    _data = grid._data;
    _cellSize = new int[_geom.getSize()[0] * _geom.getSize()[1]];
    for(int i=0; i<_cellSize.length; i++) _cellSize[i] = 1;
  }
  
  public MultiGrid(Geometry geom, int[] cellSize, float[] data) {
    _geom = geom;
    _cellSize = cellSize;
    _data = data;   
  }
  
   public void printCellSize() {
    for (int i = _geom.getSize()[1] - 1; i >= 0; i--) {
      for (int j = 0; j < _geom.getSize()[0]; j++) {
        print(" " + _cellSize[i*_geom.getSize()[0] + j] + " ");
      }
      println();
    }
    println();
  }
  
  public int getCellSize(int index) {
    return _cellSize[index];
  }
  
  public float dxx(MGIterator it) {
    if(_cellSize[it.getValue()] != 0) {
      float h = float(_cellSize[it.getValue()]);
      if(_cellSize[it.left().getValue()] == h/2) {
        if(_cellSize[it.right().getValue()] == h/2) {
          return 16/9 * (0.5 * (getCell(it.left()) + getCell(it.left().top())) - 2.0 * getCell(it) + 0.5 * (getCell(it.right()) + getCell(it.right().top()))) / (h*h);
        }
        else if(_cellSize[it.right().getValue()] == h) {
          // from matlab dxx = 2 * (4*(4*vl - 7*vc + 3*vr))/(21*h^2)
          return 8/21 * (4.0 * 0.5 * (getCell(it.left()) + getCell(it.left().top())) - 7.0 * getCell(it) + 3.0 * getCell(it.right())) / (h*h);
        }
        else {
          // from matlab dxx = 2 * (8*(2*vl - 3*vc + vr))/(27*h^2)
          return 16/27 * (2.0 * 0.5 * (getCell(it.left()) + getCell(it.left().top())) - 3.0 * getCell(it) + getCell(it.right())) / (h*h);
        }
      }
      else if(_cellSize[it.left().getValue()] == h) {
        if(_cellSize[it.right().getValue()] == h/2) {
          return 8/21 * (3.0 * getCell(it.left()) - 7.0 * getCell(it) + 4.0 * 0.5 * (getCell(it.right()) + getCell(it.right().top()))) / (h*h);
        }
        else if(_cellSize[it.right().getValue()] == h) {
          return (getCell(it.right()) - 2.0 * getCell(it) + getCell(it.left())) / h / h;
        }
        else {
          // from matlab dxx = 2 * (2*(3*vl - 5*vc + 2*vr))/(15*h^2)
          return 4/15 * (3.0 * getCell(it.left()) - 5.0 * getCell(it) + 2.0 * getCell(it.right())) / (h*h);
        }
      }
      else {
        if(_cellSize[it.right().getValue()] == h/2) {
          return 16/27 * (getCell(it.left()) - 3.0 * getCell(it) + 2.0 * 0.5 * (getCell(it.right()) + getCell(it.right().top()))) / (h*h);
        }
        else if(_cellSize[it.right().getValue()] == h) {
          return 4/15 * (2.0 * getCell(it.left()) - 5.0 * getCell(it) + 3.0 * getCell(it.right())) / (h*h);
        }
        else {
          return 4/9 * (0.5 * (getCell(it.left()) + getCell(it.left().top())) - 2.0 * getCell(it) + 0.5 * (getCell(it.right()) + getCell(it.right().top()))) / (h*h);
        }
      }
    }
    else return 0.0;
  }
  
  public float dyy(MGIterator it) {
    if(_cellSize[it.getValue()] != 0) {
      float h = float(_cellSize[it.getValue()]);
      if(_cellSize[it.top().getValue()] == h/2) {
        if(_cellSize[it.down().getValue()] == h/2) {
          return 16/9 * (0.5 * (getCell(it.top()) + getCell(it.top().top())) - 2.0 * getCell(it) + 0.5 * (getCell(it.down()) + getCell(it.down().top()))) / (h*h);
        }
        else if(_cellSize[it.down().getValue()] == h) {
          // from matlab dxx = 2 * (4*(4*vl - 7*vc + 3*vr))/(21*h^2)
          return 8/21 * (4.0 * 0.5 * (getCell(it.top()) + getCell(it.top().top())) - 7.0 * getCell(it) + 3.0 * getCell(it.down())) / (h*h);
        }
        else {
          // from matlab dxx = 2 * (8*(2*vl - 3*vc + vr))/(27*h^2)
          return 16/27 * (2.0 * 0.5 * (getCell(it.top()) + getCell(it.top().top())) - 3.0 * getCell(it) + getCell(it.down())) / (h*h);
        }
      }
      else if(_cellSize[it.top().getValue()] == h) {
        if(_cellSize[it.down().getValue()] == h/2) {
          return 8/21 * (3.0 * getCell(it.top()) - 7.0 * getCell(it) + 4.0 * 0.5 * (getCell(it.down()) + getCell(it.down().top()))) / (h*h);
        }
        else if(_cellSize[it.down().getValue()] == h) {
          return (getCell(it.down()) - 2.0 * getCell(it) + getCell(it.top())) / h / h;
        }
        else {
          // from matlab dxx = 2 * (2*(3*vl - 5*vc + 2*vr))/(15*h^2)
          return 4/15 * (3.0 * getCell(it.top()) - 5.0 * getCell(it) + 2.0 * getCell(it.down())) / (h*h);
        }
      }
      else {
        if(_cellSize[it.down().getValue()] == h/2) {
          return 16/27 * (getCell(it.top()) - 3.0 * getCell(it) + 2.0 * 0.5 * (getCell(it.down()) + getCell(it.down().top()))) / (h*h);
        }
        else if(_cellSize[it.down().getValue()] == h) {
          return 4/15 * (2.0 * getCell(it.top()) - 5.0 * getCell(it) + 3.0 * getCell(it.down())) / (h*h);
        }
        else {
          return 4/9 * (0.5 * (getCell(it.top()) + getCell(it.top().top())) - 2.0 * getCell(it) + 0.5 * (getCell(it.down()) + getCell(it.down().top()))) / (h*h);
        }
      }
    }
    else return 0.0;
  }
  
  public MultiGrid restrict(int resSize) {
    int[] newCellSize = Arrays.copyOf(_cellSize, _cellSize.length);
    float[] newData = Arrays.copyOf(_data, _data.length);
    
    MGInteriorIterator it = new MGInteriorIterator(_geom, new MultiGrid(_geom, newCellSize, newData), resSize);
    for(it.first(); it.isValid(); it.next()) {
      if(_cellSize[it.top().getValue()] == resSize && _geom.getFlag()[it.top().getValue()] == ' ' 
         && _cellSize[it.right().getValue()] == resSize & _geom.getFlag()[it.right().getValue()] == ' ' 
         && _cellSize[it.top().right().getValue()] == resSize && _geom.getFlag()[it.top().right().getValue()] == ' '
         && _cellSize[it.left().getValue()] == resSize
         && _cellSize[it.left().top().getValue()] == resSize
         && _cellSize[it.down().getValue()] == resSize
         && _cellSize[it.down().right().getValue()] == resSize
         && _cellSize[it.top().top().getValue()] == resSize
         && _cellSize[it.top().top().right().getValue()] == resSize
         && _cellSize[it.right().right().getValue()] == resSize
         && _cellSize[it.right().right().top().getValue()] == resSize) {
        
         write2Cell(it, 0.25*(getCell(it) + getCell(it.top()) + getCell(it.right()) + getCell(it.top().right())));
         write2Cell(it.getValue()+2*resSize-1, 0.25*(getCell(it) + getCell(it.top()) + getCell(it.right()) + getCell(it.top().right())));
         write2Cell(it.getValue()+2*resSize*(_geom.getSize()[0]-1), 0.25*(getCell(it) + getCell(it.top()) + getCell(it.right()) + getCell(it.top().right())));
        
         int index = it.getValue();
         int indexTop = it.top().getValue();
         int indexRight = it.right().getValue();
         int indexTopRight = it.top().right().getValue();
        
         newCellSize[index] = 2*resSize;
         newCellSize[indexTop] = 0;
         newCellSize[indexRight] = 0;
         newCellSize[indexTopRight] = 0;
      }
    }
    
    MGBoundaryIterator bit = new MGBoundaryIterator(_geom, new MultiGrid(_geom, newCellSize, newData), resSize);

    bit.setBoundary(0);
    for(bit.first(); bit.isValid(); bit.next()) {
      if(_cellSize[bit.top().getValue()] == resSize & _geom.getFlag()[bit.top().getValue()] == _geom.getFlag()[bit.getValue()]
         && _cellSize[bit.left().getValue()] == resSize
         && _cellSize[bit.left().top().getValue()] == resSize) {
         
         //write2Cell(bit, 0.5*(getCell(bit) + getCell(bit.top())));
        
         int index = bit.getValue();
         int indexRight = bit.top().getValue();
        
         newCellSize[index] = 2*resSize;
         newCellSize[indexRight] = 0;
       }
    }
    
    bit.setBoundary(1);
    for(bit.first(); bit.isValid(); bit.next()) {
      if(_cellSize[bit.right().getValue()] == resSize & _geom.getFlag()[bit.right().getValue()] == _geom.getFlag()[bit.getValue()]
         && _cellSize[bit.top().getValue()] == resSize
         && _cellSize[bit.top().right().getValue()] == resSize) {
         
         //write2Cell(bit, 0.5*(getCell(bit) + getCell(bit.right())));
         
         int index = bit.getValue();
         int indexRight = bit.right().getValue();
        
         newCellSize[index] = 2*resSize;
         newCellSize[indexRight] = 0;
       }
    }
    
    bit.setBoundary(2);
    for(bit.first(); bit.isValid(); bit.next()) {
      if(_cellSize[bit.top().getValue()] == resSize & _geom.getFlag()[bit.top().getValue()] == _geom.getFlag()[bit.getValue()]
         && _cellSize[bit.right().getValue()] == resSize
         && _cellSize[bit.right().top().getValue()] == resSize) {
         
         //write2Cell(bit, 0.5*(getCell(bit) + getCell(bit.top())));
        
         int index = bit.getValue();
         int indexRight = bit.top().getValue();
        
         newCellSize[index] = 2*resSize;
         newCellSize[indexRight] = 0;
       }
    }
    
    bit.setBoundary(3);
    for(bit.first(); bit.isValid(); bit.next()) {
      if(_cellSize[bit.right().getValue()] == resSize & _geom.getFlag()[bit.right().getValue()] == _geom.getFlag()[bit.getValue()]
         && _cellSize[bit.down().getValue()] == resSize
         && _cellSize[bit.down().right().getValue()] == resSize) {
         
         //write2Cell(bit, 0.5*(getCell(bit) + getCell(bit.right())));
        
         int index = bit.getValue();
         int indexRight = bit.right().getValue();
        
         newCellSize[index] = 2*resSize;
         newCellSize[indexRight] = 0;
       }
    }
    _geom.updateP(new MultiGrid(_geom, newCellSize, newData));
    
    return new MultiGrid(_geom, newCellSize, newData);
  }
  
  public MultiGrid interpolate(int intSize) {
    int[] newCellSize = Arrays.copyOf(_cellSize, _cellSize.length);
    float[] newData = Arrays.copyOf(_data, _data.length);
    
    MGInteriorIterator it = new MGInteriorIterator(_geom, new MultiGrid(_geom, newCellSize, newData), intSize);
    for(it.first(); it.isValid(); it.next()) {
       if(_cellSize[it.getValue()] == intSize && _geom.getFlag()[it.getValue()] == ' ') {
         
         newCellSize[it.getValue()] = intSize/2;
         newCellSize[it.getValue()+intSize/2] = intSize/2;
         newCellSize[it.getValue()+intSize/2*_geom.getSize()[0]] = intSize/2;
         newCellSize[it.getValue()+intSize/2*(1+_geom.getSize()[0])] = intSize/2;
         
         write2Cell(it.top(), getCell(it));
         write2Cell(it.right(), getCell(it));
         write2Cell(it.top().right(), getCell(it));
         
         write2Cell(it.getValue()+intSize/2-1, getCell(it));
         write2Cell(it.top().getValue()+intSize/2-1, getCell(it));
         write2Cell(it.right().getValue()+intSize/2-1, getCell(it));
         write2Cell(it.top().right().getValue()+intSize/2-1, getCell(it));
         
         write2Cell(it.getValue()+intSize/2*(_geom.getSize()[0]-1), getCell(it));
         write2Cell(it.top().getValue()+intSize/2*(_geom.getSize()[0]-1), getCell(it));
         write2Cell(it.right().getValue()+intSize/2*(_geom.getSize()[0]-1), getCell(it));
         write2Cell(it.top().right().getValue()+intSize/2*(_geom.getSize()[0]-1), getCell(it));
       }
    }
    
    MGBoundaryIterator bit = new MGBoundaryIterator(_geom, new MultiGrid(_geom, newCellSize, newData), intSize);

    bit.setBoundary(0);
    for(bit.first(); bit.isValid(); bit.next()) {
      if(_cellSize[bit.getValue()] == intSize) {
        
        newCellSize[bit.getValue()] = intSize/2;
        newCellSize[bit.getValue() + intSize/2*_geom.getSize()[0]] = intSize/2;
      }
    }
    
    bit.setBoundary(1);
    for(bit.first(); bit.isValid(); bit.next()) {
      if(_cellSize[bit.getValue()] == intSize) {
        
        newCellSize[bit.getValue()] = intSize/2;
        newCellSize[bit.getValue() + intSize/2] = intSize/2;
      }
    }
    
    bit.setBoundary(2);
    for(bit.first(); bit.isValid(); bit.next()) {
      if(_cellSize[bit.getValue()] == intSize) {
        
        newCellSize[bit.getValue()] = intSize/2;
        newCellSize[bit.getValue() + intSize/2*_geom.getSize()[0]] = intSize/2;
      }
    }
    
    bit.setBoundary(3);
    for(bit.first(); bit.isValid(); bit.next()) {
      if(_cellSize[bit.getValue()] == intSize) {
        
        newCellSize[bit.getValue()] = intSize/2;
        newCellSize[bit.getValue() + intSize/2] = intSize/2;
      }
    }
    
    return new MultiGrid(_geom, newCellSize, newData);
  }
}