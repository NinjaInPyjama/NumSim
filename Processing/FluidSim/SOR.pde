class SOR extends Solver {
  private float _omega;
  
  public SOR(Geometry geom) {
    _geom = geom;
    _omega = 2.0 / (1.0 + sin(PI*_geom.getMesh()[0]));
  }
  
  public SOR(Geometry geom, float omega) {
    _geom = geom;
    _omega = omega;
  }
  
  public float cycle(Grid grid, Grid rhs) {
    InteriorIterator it = new InteriorIterator(_geom);
    float dx = _geom.getMesh()[0];
    float dy = _geom.getMesh()[1];
    float factor = 1.0/(2.0/(dx*dx) + 2.0/(dy*dy));
    
    float total_res = 0.0;
    float local_res = 0.0;
    

    for(it.first(); it.isValid(); it.next()){
      local_res = localRes(it, grid, rhs);
      total_res += local_res*local_res;
      
      grid.add2Cell(it, _omega * local_res * factor);
    }
    
    _geom.updateP(grid);

    return total_res / (_geom.getSize()[0] * _geom.getSize()[1]);
  }
  
}