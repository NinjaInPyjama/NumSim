class MultiGridSolver extends Solver {
  int _smoothCycles;
  
  public MultiGridSolver(Geometry geom) {
    _geom = geom;
    _smoothCycles = 3;
  }
  
  public MultiGridSolver(Geometry geom, int smoothCycles) {
    _geom = geom;
    _smoothCycles = smoothCycles;
  }
  
  private float localRes(MGIterator it, MultiGrid grid, MultiGrid rhs) {
    return grid.dxx(it) + grid.dyy(it) - rhs.getCell(it);
  }
  
  private void smooth(MultiGrid grid, MultiGrid rhs) {
    
  }
  
  public void MGCycle(Grid grid, Grid rhs) {
    MGCycle(new MultiGrid(grid), new MultiGrid(rhs), 1);
  }
  
  private void MGCycle(MultiGrid grid, MultiGrid rhs, int cellSize) {
    for(int i=0; i<_smoothCycles; i++) smooth(grid, rhs); //pre-smoothing
    if(cellSize < 0.5*min(_geom.getSize()[0], _geom.getSize()[1])) {
      // TODO: solver Cycle
    }
    for(int i=0; i<_smoothCycles; i++) smooth(grid, rhs); //post-smoothing
  }
} //<>//