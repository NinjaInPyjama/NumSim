class PathLine extends ParticleLine {
  
  public PathLine(float[] pos) {
    _particles = new ArrayList<Particle>();
    _particles.add(new Particle(pos));
}

  public void timeStep(float dt, Grid u, Grid v) {
      _particles.add(0, new Particle(_particles.get(0).getPos()));
      _particles.get(0).timeStep(dt, u, v);
  }
  
  public void draw(float width, float height, Geometry geom) {
    stroke(0,0,0);
    strokeWeight(2);
    float[] pos1;
    float[] pos2;
    float[] geomLength = geom.getLength();
    for(int i=0; i<_particles.size()-1; i++) {
      pos1 = _particles.get(i).getPos();
      pos2 = _particles.get(i+1).getPos();
      line(pos1[0]/geomLength[0] * width, (1 - pos1[1]/geomLength[1]) * height, pos2[0]/geomLength[0] * width, (1 - pos2[1]/geomLength[1]) * height);
    }
  }
}