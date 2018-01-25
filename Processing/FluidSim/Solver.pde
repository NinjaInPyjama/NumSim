class Solver {
  protected Geometry _geom;
  
  public Solver(){}
  
  public Solver(Geometry geom) {
    _geom = geom;
  }
  
  public float cycle(Grid grid, Grid rhs) {
    return 0.0;
  }
  
  protected float localRes(Iterator it, Grid grid, Grid rhs) {
    return grid.dxx(it) + grid.dyy(it) - rhs.getCell(it);
  }
}