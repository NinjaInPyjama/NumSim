class Particle {
  private float[] _pos;
  
  public Particle (float[] pos) {
    _pos = Arrays.copyOf(pos,2);
  }
  
  public void timeStep(float dt, Grid u, Grid v) {
    _pos[0] = _pos[0] + dt * u.interpolate(_pos);
    _pos[1] = _pos[1] + dt * v.interpolate(_pos);
  }
  
  public float[] getPos() {
    return _pos;
  }
  
  public void draw(float width, float height, Geometry geom) {
    noStroke();
    fill(0,0,0);
    ellipse(_pos[0]/geom.getLength()[0] * width, (1 - _pos[1]/geom.getLength()[1]) * height, 3, 3);
  }
}