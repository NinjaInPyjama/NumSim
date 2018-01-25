class Visu {
  private float _width;
  private float _height;
  private int[] _res;
  
  private Geometry _geom;
  private Grid _grid;
  private float[] _data;
  
  private boolean _showPathLines;
  private boolean _showStreakLines;
  private ArrayList<PathLine> _pathLines;
  private ArrayList<StreakLine> _streakLines;
  
  private float _maxValue;
  private float _minValue;
  
  private boolean _showNiveau;
  private int _amountNiveauLines;
  
  public Visu(float width, float height, Geometry _geom) {
    _width = width;
    _height = height;
    _res = new int[]{400, 400};
    
    _data = new float[_res[0] * _res[1]];
    _geom = geom;
    _grid = new Grid(_geom);
    
    _pathLines = new ArrayList<PathLine>();
    _streakLines = new ArrayList<StreakLine>();
    
    _maxValue = 0.0;
    _minValue = 0.0;
    
    _showNiveau = false;
    _amountNiveauLines = 0;
  }
  
  public Visu(float width, float height, Grid grid) {
    _width = width;
    _height = height;
    _res = new int[]{400, 400};
    
    _data = new float[_res[0] * _res[1]];
    _grid = grid;
    _geom = grid.getGeometry();
    setGrid(_grid);
    
    _showPathLines = false;
    _showStreakLines = false;
    _pathLines = new ArrayList<PathLine>();
    _streakLines = new ArrayList<StreakLine>();
    
    _maxValue = 0.0;
    _minValue = 0.0;
    
    _showNiveau = false;
    _amountNiveauLines = 0;
  }
  
  public float getMaxValue() {
    return _maxValue;
  }
  
  public void setGrid(Grid grid) {
    _geom = grid.getGeometry();
    
    for(int j=0; j<_res[1]; j++) {
      for(int i=0; i<_res[0]; i++) {
        _data[j * _res[0] + i] = grid.interpolate(new float[]{(i+0.5) * _geom.getLength()[0] / _res[0], (j+0.5) * _geom.getLength()[1] / _res[1]});
      }
    }
    
    setMinMax();
  }
  
  public void setPathLine(ArrayList<PathLine> pathLines) {
    _pathLines = pathLines;
    _showPathLines = true;
    _showStreakLines = false;
  }
  
  public void setStreakLine(ArrayList<StreakLine> streakLines) {
    _streakLines = streakLines;
    _showPathLines = false;
    _showStreakLines = true;
  }
  
  public void showPathLines(boolean show) {
    _showPathLines = show;
  }
  
  public void showStreakLines(boolean show) {
    _showStreakLines = show;
  }
  
  public void setResolution(int[] res) {
    _res = res;
    _data = new float[_res[0] * _res[1]];
    setGrid(_grid);
  }
  
  public void draw() {
    noStroke();
    //stroke(0);
    rectMode(CORNER);
    colorMode(HSB, 360);
    for(int j=0; j<_res[1]; j++) {
      for(int i=0; i<_res[0]; i++) {
        fill(getColor(_data[j * _res[0] + i]),360,360);
        rect(i * _width / _res[0], _height - (j+1) * _height / _res[1],  _width / _res[0],  _height / _res[1]);
      }
    }
    if(_showNiveau) drawNiveau();
    if(_showPathLines) for(PathLine p : _pathLines) p.draw(_width, _height, _geom);
    if(_showStreakLines) for(StreakLine s : _streakLines) s.draw(_width, _height, _geom);
  }
  
  private void drawNiveau() {
    for(float i=_minValue; i<_maxValue; i += (_maxValue - _minValue) / _amountNiveauLines) getNiveauLine(i);
  }
  
  private void getNiveauLine(float c) {
    int pointID = 0;
    float[][] points = new float[2][2];
    float x = 0.0;
    float y = 0.0;
    int index = 0;
    
    stroke(0);
    strokeWeight(1);
    for(int j=0; j<_res[1] - 1; j++) {
      for(int i=0; i<_res[0] - 1; i++) {
        points = new float[2][2];
        pointID = 0;
        x = (i + 0.5) * _width / _res[0];
        y = _height - (j + 0.5) * _height / _res[1];
        index = j * _res[0] + i;
        
        if(pointID < 2 && _data[index] < c) {
          if(pointID < 2 && _data[index + 1] > c) {
            points[pointID] = new float[]{x + (c - _data[index]) / (_data[index + 1] - _data[index]) * _width / _res[0], y};
            pointID++;
          }
          if(pointID < 2 && _data[index + _res[0]] > c) {
            points[pointID] = new float[]{x, y - (c - _data[index]) / (_data[index + _res[0]] - _data[index]) * _height / _res[1]};
            pointID++;
          }
        }
        if(pointID < 2 && _data[index + 1] < c) {
          if(pointID < 2 && _data[index] > c) {
            points[pointID] = new float[]{x + (1 - (c - _data[index + 1]) / (_data[index] - _data[index + 1])) * _width / _res[0], y};
            pointID++;
          }
          if(pointID < 2 && _data[index + _res[0] + 1] > c) {
            points[pointID] = new float[]{x + _width / _res[0], y - (c - _data[index + 1]) / (_data[index + _res[0] + 1] - _data[index + 1]) * _height / _res[1]};
            pointID++;
          }
        }
        if(pointID < 2 && _data[index + _res[0]] < c) {
          if(pointID < 2 && _data[index + _res[0] + 1] > c) {
            points[pointID] = new float[]{x + (c - _data[index + _res[0]]) / (_data[index + _res[0]+ 1] - _data[index + _res[0]]) * _width / _res[0], y - _height / _res[1]};
            pointID++;
          }
          if(pointID < 2 && _data[index] > c) {
            points[pointID] = new float[]{x, y - (1 - (c - _data[index + _res[0]]) / (_data[index] - _data[index + _res[0]])) * _height / _res[1]};
            pointID++;
          }
        }
        if(pointID < 2 && _data[index + _res[0] + 1] < c) {
          if(pointID < 2 && _data[index + _res[0]] > c) {
            points[pointID] = new float[]{x + (1 - (c - _data[index + _res[0] + 1]) / (_data[index + _res[0]] - _data[index + _res[0] + 1])) * _width / _res[0], y - _height / _res[1]};
            pointID++;
          }
          if(pointID < 2 && _data[index + 1] > c) {
            points[pointID] = new float[]{x + _width / _res[0], y - (1 - (c - _data[index + _res[0] + 1]) / (_data[index + 1] - _data[index + _res[0] + 1])) * _height / _res[1]};
            pointID++;
          }
        }
        if(pointID == 2) line(points[0][0], points[0][1], points[1][0], points[1][1]);
      }
    }
  }
  
  public void showNiveauLines(int amount) {
    _showNiveau = true;
    _amountNiveauLines = amount;
    drawNiveau();
  }
  
  public void removeNiveauLines() {
    _showNiveau = false;
    _amountNiveauLines = 0;
  }
  
  private float getColor(float data) {
    return - 240.0 / (_maxValue - _minValue) * data + 240.0 / (_maxValue - _minValue) *_maxValue;
  }
  
  private void setMinMax() {
    _maxValue = _data[0];
    _minValue = _data[0];
    for(int i=1; i<_data.length; i++) {
      _maxValue = max(_maxValue,_data[i]);
      _minValue = min(_minValue,_data[i]);
    }
    //_minValue = -0.0065394863;
    //_maxValue = 0.42753205;
  }
}