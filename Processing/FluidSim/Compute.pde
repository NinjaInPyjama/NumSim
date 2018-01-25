class Compute {
  private float _t;
  private float _dtlimit;
  private float _epslimit;

  private Grid _u;
  private Grid _v;
  private Grid _p;

  private Grid _F;
  private Grid _G;

  // right-hand side
  private Grid _rhs;

  private SOR _solver;
  
  private ArrayList<PathLine> _pathLines;
  private ArrayList<StreakLine> _streakLines;

  private final Geometry _geom;
  private final Parameter _param;
  
  public Compute(Geometry geom, Parameter param) {
    _geom = geom;
    _param = param; 
    
    _F = new Grid(_geom, new float[]{1.0, 0.5});
    _F.initialize(0.0);
    _G = new Grid(_geom, new float[]{0.5, 1.0});
    _G.initialize(0.0);

    _u = new Grid(_geom, new float[]{1.0, 0.5});
    _u.initialize(0.0);
    _v = new Grid(_geom, new float[]{0.5, 1.0});
    _v.initialize(0.0);
    _p = new Grid(_geom, new float[]{0.5, 0.5});
    _p.initialize(0.0);

    _geom.updateU(_u);
    _geom.updateV(_v);
    _geom.updateP(_p);

    _rhs = new Grid(_geom, new float[]{0.5, 0.5});
    _rhs.initialize(0.0);

    _solver = new SOR(_geom);

    _pathLines = new ArrayList<PathLine>();
    //_pathLines.add(new PathLine(new float[]{0.01, 0.1}));
    //_pathLines.add(new PathLine(new float[]{0.01, 0.3}));
    //_pathLines.add(new PathLine(new float[]{0.01, 0.5}));
    //_pathLines.add(new PathLine(new float[]{0.01, 0.7}));
    //_pathLines.add(new PathLine(new float[]{0.01, 0.9}));
    _streakLines = new ArrayList<StreakLine>();
    //_streakLines.add(new StreakLine(new float[]{0.1, 0.1}));
    
    _t = 0.0;
    _dtlimit = _param.getDt();
    _epslimit = _param.getEps();
  }
  
  public void timeStep(boolean printInfo) {
    /*
    float dtlimit_diff = _param.getRe() / 2.0 * (_geom.getMesh()[0] * _geom.getMesh()[0]*_geom.getMesh()[1] * _geom.getMesh()[1]) 
                        / (_geom.getMesh()[0] * _geom.getMesh()[0] + _geom.getMesh()[1] * _geom.getMesh()[1]);
                        
    float dtlimit_conv_x = _dtlimit;
    float dtlimit_conv_y = _dtlimit;
    
    if(_u.getAbsMax() != 0) dtlimit_conv_x = _geom.getMesh()[0] / _u.getAbsMax();
    if(_v.getAbsMax() != 0) dtlimit_conv_y = _geom.getMesh()[1] / _v.getAbsMax();
    
    float dt = 0.9 * min(min(dtlimit_diff, _dtlimit), min(dtlimit_conv_x, dtlimit_conv_y));
    */
    float dt = _param.getDt();
    
    _geom.updateU(_u);
    _geom.updateV(_v);
    _geom.updateP(_p);
    
    
    MomentumEqu(dt);
    
    _geom.updateU(_F);
    _geom.updateV(_G);
    
    // compute rhs
    RHS(dt);

    // solver iterations
    int it = 0;
    float res = 0.0;    
    
    do {
        it++;
        res = _solver.cycle(_p, _rhs);
    } while(it<_param.getItermax() && res > _epslimit*_epslimit);
    if(printInfo) println("Solver stopped at iteration " + it + " with residual: " + sqrt(res));
  
    newVelocities(dt);    
    
    _geom.updateU(_u);
    _geom.updateV(_v);
    _geom.updateP(_p);
  
    
    for(PathLine p : _pathLines) p.timeStep(dt, _u, _v);
    for(StreakLine s : _streakLines) s.timeStep(dt, _u, _v);
    
    // save timestep
    _t += dt;
  }

  public float getTime() {
    return _t;
  }

  public Grid getU() {
    return _u.copy();
  }
  
  public Grid getV() {
    return _v.copy();
  }
  
  public Grid getP() {
    return _p.copy();
  }
  
  public Grid getRHS() {
    return _rhs.copy();
  }
  
  public ArrayList<PathLine> getPathLines() {
    return _pathLines;
  }
  
  public ArrayList<StreakLine> getStreakLines() {
    return _streakLines;
  }

  public Grid getVelocity() {
    InteriorIterator iit = new InteriorIterator(_geom);
    Grid absVel = new Grid(_geom, new float[]{0.5, 0.5});
    absVel.initialize(0.0);
    
    float u_ip = 0.0;
    float v_ip = 0.0;

    for(iit.first(); iit.isValid(); iit.next()){
      // Interpolating the velocities to center of cells
      v_ip = 0.5 * (_v.getCell(iit.down()) + _v.getCell(iit));
      u_ip = 0.5 * (_u.getCell(iit.left()) + _u.getCell(iit));
      absVel.write2Cell(iit, sqrt(v_ip*v_ip + u_ip*u_ip));
    }
      
    BoundaryIterator bit = new BoundaryIterator(_geom);

    for (bit.first(); bit.isValid(); bit.next()) {
      absVel.write2Cell(bit, 0.0);
      switch (_geom.getFlag()[bit.getValue()]) {
      case '#': // NOSLIP
        if (_geom.getFlag()[bit.top().getValue()] == ' ') absVel.write2Cell(bit, -absVel.getCell(bit.top()));
        if (_geom.getFlag()[bit.down().getValue()] == ' ') absVel.write2Cell(bit, -absVel.getCell(bit.down()));
        if (_geom.getFlag()[bit.left().getValue()] == ' ') absVel.write2Cell(bit, -absVel.getCell(bit.left())); 
        if (_geom.getFlag()[bit.right().getValue()] == ' ') absVel.write2Cell(bit, -absVel.getCell(bit.right()));
        break;
      case '-': // Horizontal SLIP
        if (_geom.getFlag()[bit.top().getValue()] == '-' || _geom.getFlag()[bit.down().getValue()] == '-') {
          if (_geom.getFlag()[bit.right().getValue()] == ' ') absVel.write2Cell(bit, absVel.getCell(bit.right()));
          else absVel.write2Cell(bit, absVel.getCell(bit.left()));
        }
        else {
          if (_geom.getFlag()[bit.top().getValue()] == ' ') absVel.write2Cell(bit, -absVel.getCell(bit.top()));
          if (_geom.getFlag()[bit.down().getValue()] == ' ') absVel.write2Cell(bit, -absVel.getCell(bit.down()));
          if (_geom.getFlag()[bit.left().getValue()] == ' ') absVel.write2Cell(bit, -absVel.getCell(bit.left()));
          if (_geom.getFlag()[bit.right().getValue()] == ' ') absVel.write2Cell(bit, -absVel.getCell(bit.right()));
        }
        break;
      case '|': // Vertical SLIP
        if (_geom.getFlag()[bit.left().getValue()] == '|' || _geom.getFlag()[bit.right().getValue()] == '|') {
          if (_geom.getFlag()[bit.top().getValue()] == ' ') absVel.write2Cell(bit, absVel.getCell(bit.top()));
          else absVel.write2Cell(bit, absVel.getCell(bit.down()));
        }
        else {
          if (_geom.getFlag()[bit.top().getValue()] == ' ') absVel.write2Cell(bit, -absVel.getCell(bit.top()));
          if (_geom.getFlag()[bit.down().getValue()] == ' ') absVel.write2Cell(bit, -absVel.getCell(bit.down()));
          if (_geom.getFlag()[bit.left().getValue()] == ' ') absVel.write2Cell(bit, -absVel.getCell(bit.left()));
          if (_geom.getFlag()[bit.right().getValue()] == ' ')absVel.write2Cell(bit, -absVel.getCell(bit.right()));
        }
        break;
      case 'O': // OUTFLOW
        if (_geom.getFlag()[bit.top().getValue()] == ' ') absVel.write2Cell(bit, absVel.getCell(bit.top()));
        else if (_geom.getFlag()[bit.down().getValue()] == ' ') absVel.write2Cell(bit, absVel.getCell(bit.down()));
        else if (_geom.getFlag()[bit.left().getValue()] == ' ') absVel.write2Cell(bit,absVel.getCell(bit.left()));
        else if (_geom.getFlag()[bit.right().getValue()] == ' ')absVel.write2Cell(bit, absVel.getCell(bit.right()));
        break;
      case 'V': // Vertical INFLOW
        if (_geom.getFlag()[bit.top().getValue()] == ' ') absVel.write2Cell(bit, 0.5 * (_u.getCell(bit.left()) + _u.getCell(bit)));
        else if (_geom.getFlag()[bit.down().getValue()] == ' ') absVel.write2Cell(bit, 0.5 * (_u.getCell(bit.left()) + _u.getCell(bit)));
        else if (_geom.getFlag()[bit.left().getValue()] == ' ') absVel.write2Cell(bit, 2.0 * _u.getCell(bit.left()) - absVel.getCell(bit.left()));
        else if (_geom.getFlag()[bit.right().getValue()] == ' ') absVel.write2Cell(bit, 2.0 * _u.getCell(bit) - absVel.getCell(bit.right()));
        break;
      case 'H': // Horizontal INFLOW
        if (_geom.getFlag()[bit.top().getValue()] == ' ') absVel.write2Cell(bit, 2.0 * _v.getCell(bit) - absVel.getCell(bit.top()));
        else if (_geom.getFlag()[bit.down().getValue()] == ' ') absVel.write2Cell(bit, 2.0 * _v.getCell(bit.down()) - absVel.getCell(bit.down()));
        else if (_geom.getFlag()[bit.left().getValue()] == ' ') absVel.write2Cell(bit, 0.5 * (_v.getCell(bit.top()) + _v.getCell(bit)));
        else if (_geom.getFlag()[bit.right().getValue()] == ' ') absVel.write2Cell(bit, 0.5 * (_v.getCell(bit.top()) + _v.getCell(bit)));
        break;
      default:
        break;
      }
    }
  
    return absVel;
  }
  
  
  public Grid getVorticity() {
    Iterator it = new Iterator(_geom);
    Grid vort = new Grid(_geom, new float[]{1.0, 1.0});
    
    vort.initialize(0.0);
    
    for(it.first(); it.isValid(); it.next()){ 
        // Calculating vorticity by du/dy - dv/dx
    if( _geom.getFlag()[it.getValue()] == ' ' ) {
            vort.write2Cell(it, _u.dy_r(it) - _v.dx_r(it));
        }
    }
    return vort;
  }
  
  public Grid getStream() {
    Grid psi = new Grid(_geom, new float[]{1.0, 1.0});
    psi.initialize(0.0);
    /*
    BoundaryIterator bit = new BoundaryIterator(_geom);
    bit.setBoundary(3);
    bit.first();
    bit.next();
    for (; bit.isValid(); bit.next()) {
      psi.write2Cell(bit, _geom.getMesh()[1] * _u.getCell(bit) + psi.getCell(bit.down()));
    }
    
    bit.setBoundary(2);
    bit.first();
    bit.next();
    for (; bit.isValid(); bit.next()) {
        psi.write2Cell(bit, - _geom.getMesh()[0] * _v.getCell(bit) + psi.getCell(bit.left()));
    }
    
    InteriorIterator iit = new InteriorIterator(_geom);
    
    for(iit.first(); iit.isValid(); iit.next()){
        psi.write2Cell(iit, - _geom.getMesh()[0] * _v.getCell(iit) + psi.getCell(iit.left()));
    }
    */
    return psi;
  }
  
  private void newVelocities(float dt) {
    InteriorIterator iit = new InteriorIterator(_geom);
    
    for(iit.first(); iit.isValid(); iit.next()) {
      _u.write2Cell(iit, _F.getCell(iit) - dt*(_p.dx_r(iit)));
      _v.write2Cell(iit, _G.getCell(iit) - dt*(_p.dy_r(iit)));
    }
  }
  
  private void MomentumEqu(float dt) {
    InteriorIterator iit = new InteriorIterator(_geom);

    for(iit.first(); iit.isValid(); iit.next()) {
      float A = (_u.dxx(iit) + _u.dyy(iit)) / _param.getRe() - _u.DC_du2_x(iit, _param.getAlpha()) - _u.DC_duv_y(iit, _param.getAlpha(), _v);
      _F.write2Cell(iit, _u.getCell(iit) + dt*A);
      
      float B = (_v.dxx(iit) + _v.dyy(iit)) / _param.getRe() - _v.DC_dv2_y(iit, _param.getAlpha()) - _v.DC_duv_x(iit, _param.getAlpha(), _u);
      _G.write2Cell(iit, _v.getCell(iit) + dt*B);
    }
  }
  
  private void RHS(float dt) {
    InteriorIterator iit = new InteriorIterator(_geom);

    for(iit.first(); iit.isValid(); iit.next()) _rhs.write2Cell(iit, (_F.dx_l(iit) + _G.dy_l(iit))/dt);
  }
}