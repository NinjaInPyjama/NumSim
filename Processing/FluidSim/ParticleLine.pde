class ParticleLine {
  protected ArrayList<Particle> _particles;
  
  public void timeStep(float dt, Grid u, Grid v){}
  
  public ArrayList<Particle> getParticles() {
    return _particles;
  }
  
  public void draw() {}
}