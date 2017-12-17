#include "compute.hpp"

using namespace std;

/// Creates a compute instance with given geometry and parameter
Compute::Compute(const Geometry * geom, const Parameter * param) {
    
    _geom = geom;
    _param = param;
    
    _tmp = new Grid(_geom);
    
    _F = new Grid(_geom, multi_real_t(1.0, 0.5));
	_F->Initialize(0.0);
    _G = new Grid(_geom, multi_real_t(0.5, 1.0));
	_G->Initialize(0.0);

	_u = new Grid(_geom, multi_real_t(1.0, 0.5));
	_u->Initialize(0.0);
	_v = new Grid(_geom, multi_real_t(0.5, 1.0));
	_v->Initialize(0.0);
    _p = new Grid(_geom, multi_real_t(0.5, 0.5));
    _p->Initialize(_geom->Pressure());

	_geom->Update_U(_u);
	_geom->Update_V(_v);
	_geom->Update_P(_p);

	_rhs = new Grid(_geom, multi_real_t(0.5, 0.5));
	_rhs->Initialize(0.0);

	_solver = new SOR(_geom);
        
    _t = 0.0;
    _dtlimit = _param->Dt();
    _epslimit = _param->Eps();
    _pathline = new PathLine(multi_real_t( 0.2, 0.2));
    _streakline = new StreakLine(multi_real_t( 0.2, 0.8));
}

/// Deletes all grids
Compute::~Compute() {}


/// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual
// etc.)
void Compute::TimeStep(bool printInfo) {
    // see script p. 23
    
    real_t dtlimit_diff = _param->Re()/2.0 * (_geom->Mesh()[0]*_geom->Mesh()[0]*_geom->Mesh()[1]*_geom->Mesh()[1])/(_geom->Mesh()[0]*_geom->Mesh()[0]+_geom->Mesh()[1]*_geom->Mesh()[1]);
    
    // stability  condition induced by the convection operator
    real_t dtlimit_conv_x = _dtlimit;
    real_t dtlimit_conv_y = _dtlimit;
    real_t AbsMax_v = _v->AbsMax();
    real_t AbsMax_u = _u->AbsMax();
    
    if(AbsMax_u != 0)
        dtlimit_conv_x = _geom->Mesh()[0]/AbsMax_u;
    
    if(AbsMax_v != 0)
        dtlimit_conv_y = _geom->Mesh()[1]/AbsMax_v;
    
    
    
    // minimum of all time limits
    real_t dt = 0.2*std::min(std::min(dtlimit_diff, _dtlimit),std::min(dtlimit_conv_x, dtlimit_conv_y));
    
    //dt = 0.005;
	// update boundary values
	_geom->Update_U(_u);
	_geom->Update_V(_v);
	_geom->Update_P(_p);

    // compute 'preliminary' velocities and setting boundary values
    MomentumEqu(dt);
	_geom->Update_U(_F);
	_geom->Update_V(_G);

    // compute rhs
    RHS(dt);

    // solver iterations
    index_t itermax = _param->IterMax();
    index_t it = 0;
    real_t res = 0.0;
    do {
        it++;
        res = _solver->Cycle(_p, _rhs);
		//if(printInfo) std::cout << "Residual at iteration " << it << ": " << res << std::endl;
    } while(it<itermax && res>_epslimit);
    if(printInfo) std::cout << "Solver stopped at iteration " << it << " with residual: " << res << " with dt conv x: " << dtlimit_conv_x << " with dt diff: " << dtlimit_diff << " t: " << _t << std::endl;
	
    // compute new velocities
    NewVelocities(dt);

	// udating boundary values (to be consistent when saving vtks)
	_geom->Update_U(_u);
	_geom->Update_V(_v);
	_geom->Update_P(_p);
    
    _pathline->TimeStep(dt, _u, _v);
    _streakline->TimeStep(dt, _u, _v);
    
    //_pathline->print();
	//_streakline->print();
	// save timestep
    _t += dt;
}

/// Returns the simulated time in total
const real_t & Compute::GetTime() const {
	return _t;
}

/// Returns the simulated PathLine
const PathLine * Compute::GetPathLine() const {
	return _pathline;
}

/// Returns the simulated StreakLines
const StreakLine * Compute::GetStreakLine() const {
	return _streakline;
}


/// Returns the pointer to U
const Grid * Compute::GetU() const {
	return _u;
}

/// Returns the pointer to V
const Grid * Compute::GetV() const {
	return _v;
}

/// Returns the pointer to P
const Grid * Compute::GetP() const {
	return _p;
}

/// Returns the pointer to RHS
const Grid * Compute::GetRHS() const {
	return _rhs;
}


/// Computes and returns the absolute velocity
const Grid * Compute::GetVelocity() {
    InteriorIterator iit = InteriorIterator(_geom);
    Grid * abs_vel = new Grid(_geom,multi_real_t(0.5,0.5));
	real_t u_ip = 0.0; // storage for interpolated u to center of cells
	real_t v_ip = 0.0; // storage for interpolated v to center of cells


    for(iit.First(); iit.Valid(); iit.Next()){
        // Interpolating the velocities to center of cells
		v_ip = (_v->Cell(iit.Down()) + _v->Cell(iit)) / 2.0;
		u_ip = (_u->Cell(iit.Left()) + _u->Cell(iit)) / 2.0;
        abs_vel->Cell(iit) = sqrt(v_ip*v_ip + u_ip*u_ip);
    }

	BoundaryIterator bit = BoundaryIterator(_geom);

	for (bit.First(); bit.Valid(); bit.Next()) {
		abs_vel->Cell(bit) = 0.0;
		switch (_geom->Flag()[bit.Value()]) {
		case '#': // NOSLIP
			if (_geom->Flag()[bit.Top().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Top());
			if (_geom->Flag()[bit.Down().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Down());
			if (_geom->Flag()[bit.Left().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Left()); 
			if (_geom->Flag()[bit.Right().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Right());
			break;
		case '-': // Horizontal SLIP
			if (_geom->Flag()[bit.Top().Value()] == '-' || _geom->Flag()[bit.Down().Value()] == '-') {
				if (_geom->Flag()[bit.Right().Value()] == ' ') abs_vel->Cell(bit) = abs_vel->Cell(bit.Right());
				else abs_vel->Cell(bit) = abs_vel->Cell(bit.Left());
			}
			else {
				if (_geom->Flag()[bit.Top().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Top());
				if (_geom->Flag()[bit.Down().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Down());
				if (_geom->Flag()[bit.Left().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Left());
				if (_geom->Flag()[bit.Right().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Right());
			}
			break;
		case '|': // Vertical SLIP
			if (_geom->Flag()[bit.Left().Value()] == '|' || _geom->Flag()[bit.Right().Value()] == '|') {
				if (_geom->Flag()[bit.Top().Value()] == ' ') abs_vel->Cell(bit) = abs_vel->Cell(bit.Top());
				else abs_vel->Cell(bit) = abs_vel->Cell(bit.Down());
			}
			else {
				if (_geom->Flag()[bit.Top().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Top());
				if (_geom->Flag()[bit.Down().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Down());
				if (_geom->Flag()[bit.Left().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Left());
				if (_geom->Flag()[bit.Right().Value()] == ' ') abs_vel->Cell(bit) = -abs_vel->Cell(bit.Right());
			}
			break;
		case 'O': // OUTFLOW
			if (_geom->Flag()[bit.Top().Value()] == ' ') abs_vel->Cell(bit) = abs_vel->Cell(bit.Top());
			else if (_geom->Flag()[bit.Down().Value()] == ' ') abs_vel->Cell(bit) = abs_vel->Cell(bit.Down());
			else if (_geom->Flag()[bit.Left().Value()] == ' ') abs_vel->Cell(bit) = abs_vel->Cell(bit.Left());
			else if (_geom->Flag()[bit.Right().Value()] == ' ') abs_vel->Cell(bit) = abs_vel->Cell(bit.Right());
			break;
		case 'V': // Vertical INFLOW
			if (_geom->Flag()[bit.Top().Value()] == ' ') abs_vel->Cell(bit) = (_u->Cell(bit.Left()) + _u->Cell(bit)) / 2.0;
			else if (_geom->Flag()[bit.Down().Value()] == ' ') abs_vel->Cell(bit) = (_u->Cell(bit.Left()) + _u->Cell(bit)) / 2.0;
			else if (_geom->Flag()[bit.Left().Value()] == ' ') abs_vel->Cell(bit) = 2.0 * _u->Cell(bit.Left()) - abs_vel->Cell(bit.Left());
			else if (_geom->Flag()[bit.Right().Value()] == ' ') abs_vel->Cell(bit) = 2.0 * _u->Cell(bit) - abs_vel->Cell(bit.Right());
			break;
		case 'H': // Horizontal INFLOW
			if (_geom->Flag()[bit.Top().Value()] == ' ') abs_vel->Cell(bit) = 2.0 * _v->Cell(bit) - abs_vel->Cell(bit.Top());
			else if (_geom->Flag()[bit.Down().Value()] == ' ') abs_vel->Cell(bit) = 2.0 * _v->Cell(bit.Down()) - abs_vel->Cell(bit.Down());
			else if (_geom->Flag()[bit.Left().Value()] == ' ') abs_vel->Cell(bit) = (_v->Cell(bit.Top()) + _v->Cell(bit)) / 2.0;
			else if (_geom->Flag()[bit.Right().Value()] == ' ') abs_vel->Cell(bit) = (_v->Cell(bit.Top()) + _v->Cell(bit)) / 2.0;
			break;
		default:
			break;
		}
	}
    return abs_vel;
}

/// Computes and returns the vorticity
const Grid * Compute::GetVorticity() {
    InteriorIterator iit = InteriorIterator(_geom);
    Grid * vort = new Grid(_geom);
	vort->Initialize(0.0);
    multi_real_t cell_center = multi_real_t(0.5, 0.5);

    // Creating grids of derivatives of u (in y-dim) and v (in x-dim) (for interpolation reasons)
    Grid * du_dy = new Grid(_geom);
    Grid * dv_dx = new Grid(_geom);
    for(iit.First();iit.Valid();iit.Next()){
        du_dy->Cell(iit) = _u->dy_central(iit);
        dv_dx->Cell(iit) = _v->dx_central(iit);
    }
    
    for(iit.First();iit.Valid();iit.Next()){ 
        // Calculating vorticity by dv/dx - du/dy
        multi_index_t cell_pos = iit.Pos();
		vort->Cell(iit) = 0.0; //  dv_dx->Interpolate(multi_real_t((real_t)cell_pos[0] + cell_center[0], (real_t)cell_pos[1] + cell_center[1]))
                                - du_dy->Interpolate(multi_real_t((real_t)cell_pos[0] + cell_center[0], (real_t)cell_pos[1]  + cell_center[1]));
    }
    return vort;
}

/// Computes and returns the stream line values
const Grid * Compute::GetStream() {
    Grid * psi = new Grid(_geom, multi_real_t(1.0, 1.0));
    psi->Initialize(0.0);
    
    int * streamset = new int[_geom->Size()[0]*_geom->Size()[1]];
    for (int i = 0;i<_geom->Size()[0]*_geom->Size()[1];i++){
        streamset[i] = 0;
    }
    bool firstvalue = false;
    Iterator it = Iterator(_geom);
    BackwardsIterator backit = BackwardsIterator(_geom);
    do {
        
    for(it.First(); it.Valid(); it.Next()){
        if(_geom->Flag()[it.Value()] == ' ' || _geom->Flag()[it.Right().Top().Value()] == ' ' || _geom->Flag()[it.Right().Value()] == ' ' || _geom->Flag()[it.Top().Value()] == ' ') {
            if(firstvalue==false){
                firstvalue = true;
                streamset[it.Value()]=1;
            }
            else if(streamset[it.Left().Value()]==1 && streamset[it.Value()]!=1){
                psi->Cell(it) = - _geom->Mesh()[0]*_v->Cell(it) + psi->Cell(it.Left());
                streamset[it.Value()]=1;
            }
            else if(streamset[it.Down().Value()]==1 && streamset[it.Value()]!=1){
                psi->Cell(it) = _geom->Mesh()[1]*_u->Cell(it) + psi->Cell(it.Down());
                streamset[it.Value()]=1;
            }
            else if(streamset[it.Right().Value()]==1 && streamset[it.Value()]!=1) {
                psi->Cell(it) = _geom->Mesh()[0]*_v->Cell(it.Right()) + psi->Cell(it.Right()) ;
                streamset[it.Value()]=1;
            }
            else if(streamset[it.Top().Value()]==1 && streamset[it.Value()]!=1){
                psi->Cell(it) = - _geom->Mesh()[1]*_u->Cell(it.Top()) + psi->Cell(it.Top()) ;
                streamset[it.Value()]=1;
            }
            else if (streamset[it.Value()]==0){
                streamset[it.Value()]=-1;
            }
        }
    }
    for(backit.First(); backit.Valid(); backit.Next()){
        if(_geom->Flag()[backit.Value()] == ' ' || _geom->Flag()[backit.Right().Top().Value()] == ' ' || _geom->Flag()[backit.Right().Value()] == ' ' || _geom->Flag()[backit.Top().Value()] == ' ') {
            if(streamset[backit.Right().Value()]==1 && streamset[backit.Value()]!=1) {
                psi->Cell(backit) = _geom->Mesh()[0]*_v->Cell(backit.Right()) + psi->Cell(backit.Right()) ;
                streamset[backit.Value()]=1;
            }
            else if(streamset[backit.Down().Value()]==1 && streamset[backit.Value()]!=1){
                psi->Cell(backit) = _geom->Mesh()[1]*_u->Cell(backit) + psi->Cell(backit.Down());
                streamset[backit.Value()]=1;
            }
            else if(streamset[backit.Left().Value()]==1 && streamset[backit.Value()]!=1){
                psi->Cell(backit) = - _geom->Mesh()[0]*_v->Cell(backit) + psi->Cell(backit.Left());
                streamset[backit.Value()]=1;
            }
            else if(streamset[backit.Top().Value()]==1 && streamset[backit.Value()]!=1){
                psi->Cell(backit) = - _geom->Mesh()[1]*_u->Cell(backit.Top()) + psi->Cell(backit.Top()) ;
                streamset[backit.Value()]=1;
            }
        }
    } std::cout << nonsetstreamset(streamset) << std::endl;
    } while (nonsetstreamset(streamset));
    
    return psi;
}

bool Compute::nonsetstreamset(const int * streamset) {
    Iterator it = Iterator(_geom);
    for(it.First(); it.Valid(); it.Next()){
        if (streamset[it.Value()]==-1) return true;
    }
    return false;
}


/// Compute the new velocites u,v
void Compute::NewVelocities(const real_t & dt) {
    InteriorIterator iit = InteriorIterator(_geom);
    
    for(iit.First();iit.Valid();iit.Next()){
		// see script, p. 18, formular (3.1)
        _u->Cell(iit) = _F->Cell(iit) - dt*(_p->dx_r(iit));
        _v->Cell(iit) = _G->Cell(iit) - dt*(_p->dy_r(iit));
    }
}

/// Compute the temporary velocites F,G
void Compute::MomentumEqu(const real_t & dt) {
    InteriorIterator iit = InteriorIterator(_geom);

    for(iit.First();iit.Valid();iit.Next()){
		// see script, p. 18, formular (3.2)
		real_t A = (_u->dxx(iit) + _u->dyy(iit)) / _param->Re() - _u->DC_du2_x(iit, _param->Alpha()) - _u->DC_duv_y(iit, _param->Alpha(), _v);
		real_t B = (_v->dxx(iit) + _v->dyy(iit)) / _param->Re() - _v->DC_dv2_y(iit, _param->Alpha()) - _v->DC_duv_x(iit, _param->Alpha(), _u);
		_F->Cell(iit) = _u->Cell(iit) + dt*A;
		_G->Cell(iit) = _v->Cell(iit) + dt*B;
    }
}

/// Compute the RHS of the poisson equation
void Compute::RHS(const real_t & dt) {
    InteriorIterator iit = InteriorIterator(_geom);

    for(iit.First();iit.Valid();iit.Next()){
		// see script, p. 19, formular (3.3)
        _rhs->Cell(iit) = (_F->dx_l(iit) + _G->dy_l(iit))/dt;
    }
}
