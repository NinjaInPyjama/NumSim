class StreakLine extends ParticleLine {
  private float[] _origin;
  
  public StreakLine(float[] pos) {
    _particles = new ArrayList<Particle>();
    _particles.add(new Particle(pos));
    _origin = pos;
  }
  
  public void timeStep(float dt, Grid u, Grid v) {
      for(Particle p : _particles) p.timeStep(dt,u,v);
      _particles.add(0, new Particle(_origin));
  }
  
  public void draw(float width, float height, Geometry geom) {
    for(Particle p : _particles) p.draw(width, height, geom);
  }
}