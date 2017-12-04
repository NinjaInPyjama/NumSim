#include "compute.hpp"
#include <mpi.h>

using namespace std;

/// Creates a compute instance with given geometry and parameter
Compute::Compute(const Geometry * geom, const Parameter * param, const Communicator *comm) {

    _geom = geom;
    _param = param;
    _comm = comm;
    
    _zeit_comm = new Zeitgeist();
    _zeit_dt = new Zeitgeist();
    _zeit_comp = new Zeitgeist();
    _zeit_res = new Zeitgeist();
    
    
    
    _tmp = new Grid(_geom);
    
    _F = new Grid(_geom, multi_real_t(1.0, 0.5));
	_F->Initialize(0.0);
    _G = new Grid(_geom, multi_real_t(0.5, 1.0));
	_G->Initialize(0.0);

	_u = new Grid(_geom, multi_real_t(1.0, 0.5));
	_u->Initialize(_geom->Velocity()[0]);
	_v = new Grid(_geom, multi_real_t(0.5, 1.0));
	_v->Initialize(_geom->Velocity()[1]);
    _p = new Grid(_geom, multi_real_t(0.5, 0.5));
    _p->Initialize(_geom->Pressure());

	_geom->Update_U(_u);
	_geom->Update_V(_v);
	_geom->Update_P(_p);

	_rhs = new Grid(_geom, multi_real_t(0.5, 0.5));
	_rhs->Initialize(0.0);

	_solver = new RedOrBlackSOR(_geom,_param->Omega());
        
    _t = 0.0;
    _dtlimit = _param->Dt();
    _epslimit = _param->Eps();
}

/// Deletes all grids
Compute::~Compute() {}


/// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual
// etc.)
void Compute::TimeStep(bool printInfo) {
    // see script p. 23
    // compute dt
    // stability  condition induced by the diffusion operator
    
    //_zeit_dt->Tic();
    
    real_t dtlimit_diff = _param->Re()/2.0 * (_geom->Mesh()[0]*_geom->Mesh()[0]*_geom->Mesh()[1]*_geom->Mesh()[1])/(_geom->Mesh()[0]*_geom->Mesh()[0]+_geom->Mesh()[1]*_geom->Mesh()[1]);
    // stability  condition induced by the convection operator
    real_t dtlimit_conv_x = _dtlimit;
    real_t dtlimit_conv_y = _dtlimit;
    real_t AbsMax_v = _comm->gatherMax(_v->AbsMax());
    real_t AbsMax_u = _comm->gatherMax(_u->AbsMax());
    
    if(AbsMax_u != 0)
        dtlimit_conv_x = _geom->Mesh()[0]/AbsMax_u;
    
    if(AbsMax_v != 0)
        dtlimit_conv_y = _geom->Mesh()[1]/AbsMax_v;
    
    // minimum of all time limits
    real_t dt = 0.9*std::min(std::min(dtlimit_diff, _dtlimit),std::min(dtlimit_conv_x, dtlimit_conv_y));
    
    
    //_zeit_dt->Tac();
    
    
    //if(_comm->getRank()==0) std::cout << "Zeitberechnung: " << _zeit->Toc() << std::endl;
    
    //if(printInfo && _comm->getRank()==0) std::cout << "Timestep dt: " << dt << std::endl;
	
    
    
    
	// update boundary values
	_geom->Update_U(_u);
	_geom->Update_V(_v);
	_geom->Update_P(_p);
    
    
    //GetStream()->print();

    // compute 'preliminary' velocities and setting boundary values
    MomentumEqu(dt);
    
	_comm->copyBoundary(_F);
	_comm->copyBoundary(_G);
    
    _geom->Update_U(_F);
	_geom->Update_V(_G);
    
    // compute rhs
    RHS(dt);

    // solver iterations
    index_t itermax = _param->IterMax();
    index_t it = 0;
    real_t res = 0.0;
    
    _comm->copyBoundary(_p);
    
    real_t zahl = 0.0;
    
    
    do {
        it++;
        //_zeit_comp->Tic();
         zahl = _solver->RedCycle(_p, _rhs);
        //_zeit_comp->Tac();
        //_zeit_comm->Tic();
        res = _comm->gatherSum(zahl);
        _comm->copyBoundary(_p);
        //_zeit_comm->Tac();
        //_zeit_comp->Tic();
        zahl = _solver->BlackCycle(_p, _rhs);
        //_zeit_comp->Tac();
        //_zeit_comm->Tic();
        res += _comm->gatherSum(zahl);
        _comm->copyBoundary(_p);
        //_zeit_comm->Tac();
        
        printTimes();
        
    } while(false && it<itermax && res>_epslimit*_epslimit);
    if(printInfo && _comm->getRank()==0) std::cout << "Solver stopped at iteration " << it << " with residual**2: " << res << std::endl;
	
    //_solver->printTimes();
    // compute new velocities
    
    NewVelocities(dt);
        
	// udating boundary values (to be consistent when saving vtks)
    
    _comm->copyBoundary(_u);
    _comm->copyBoundary(_v);
    _comm->copyBoundary(_p);
    
    
	_geom->Update_V(_v);
    _geom->Update_U(_u);
	_geom->Update_P(_p);
	
	// save timestep
    _t += dt;    
}

/// Returns the simulated time in total
const real_t & Compute::GetTime() const {
	return _t;
}

void  Compute::printTimes()  {
	std::cout << "Computation: " << _zeit_comp->Toc() << "Communication: " << _zeit_comm->Toc() << "Timestep: " << _zeit_dt->Toc() << "Residual: " << _zeit_res->Toc() << std::endl;
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
  InteriorIterator iit(_geom);
	BoundaryIterator bit(_geom);
  Grid * abs_vel = new Grid(_geom, multi_real_t(0.5, 0.5));
	
	real_t u_ip = 0.0; // storage for interpolated u to center of cells
	real_t v_ip = 0.0; // storage for interpolated v to center of cells
	
	abs_vel->Initialize(0.0);
    
  for(iit.First(); iit.Valid(); iit.Next()){
      // Interpolating the velocities to center of cells  
  v_ip = (_v->Cell(iit.Down()) + _v->Cell(iit)) / 2.0;
  u_ip = (_u->Cell(iit.Left()) + _u->Cell(iit)) / 2.0;
      abs_vel->Cell(iit) = sqrt(v_ip*v_ip + u_ip*u_ip);
  }

  bit.SetBoundary(bit.boundaryTop);
  if(_comm->isTop()) {
    for (bit.First(); bit.Valid(); bit.Next()) {
      abs_vel->Cell(bit) = 2.0 - abs_vel->Cell(bit.Down());
    }
  }
  else{
    for (bit.First(); bit.Valid(); bit.Next()) {
      v_ip = (_v->Cell(bit.Down()) + _v->Cell(bit)) / 2.0;
      u_ip = (_u->Cell(bit.Left()) + _u->Cell(bit)) / 2.0;
      abs_vel->Cell(bit) = sqrt(v_ip*v_ip + u_ip*u_ip);
    }
  }

  bit.SetBoundary(bit.boundaryBottom);
  if(_comm->isBottom()) {
    for (bit.First(); bit.Valid(); bit.Next()) {
      abs_vel->Cell(bit) = - abs_vel->Cell(bit.Top());
    }
  }
  else {
    for (bit.First(); bit.Valid(); bit.Next()) {
      v_ip = (_v->Cell(bit.Down()) + _v->Cell(bit)) / 2.0;
      u_ip = (_u->Cell(bit.Left()) + _u->Cell(bit)) / 2.0;
      abs_vel->Cell(bit) = sqrt(v_ip*v_ip + u_ip*u_ip);
    }
  }

  bit.SetBoundary(bit.boundaryRight);
  if(_comm->isRight()) {
    for (bit.First(); bit.Valid(); bit.Next()) {
      abs_vel->Cell(bit) = - abs_vel->Cell(bit.Left());
    }
  }
  else {
    for (bit.First(); bit.Valid(); bit.Next()) {
      v_ip = (_v->Cell(bit.Down()) + _v->Cell(bit)) / 2.0;
      u_ip = (_u->Cell(bit.Left()) + _u->Cell(bit)) / 2.0;
      abs_vel->Cell(bit) = sqrt(v_ip*v_ip + u_ip*u_ip);
    }
  }


  bit.SetBoundary(bit.boundaryLeft);
  if(_comm->isLeft()) {
    for (bit.First(); bit.Valid(); bit.Next()) {
      abs_vel->Cell(bit) = - abs_vel->Cell(bit.Right());
    }
  }
  else {
    for (bit.First(); bit.Valid(); bit.Next()) {
      v_ip = (_v->Cell(bit.Down()) + _v->Cell(bit)) / 2.0;
      u_ip = (_u->Cell(bit.Left()) + _u->Cell(bit)) / 2.0;
      abs_vel->Cell(bit) = sqrt(v_ip*v_ip + u_ip*u_ip);
    }
  }

  return abs_vel;
}

/// Computes and returns the vorticity
const Grid * Compute::GetVorticity() {
    InteriorIterator iit(_geom);
	BoundaryIterator bit = BoundaryIterator(_geom);
    Grid * vort = new Grid(_geom, multi_real_t(1.0, 1.0));
    
    vort->Initialize(0.0);
    
	bit.SetBoundary(bit.boundaryBottom);
	for(bit.First(); bit.Valid(); bit.Next()) {
		vort->Cell(bit) = _u->dy_r(bit) - _v->dx_r(bit);
	}

	bit.SetBoundary(bit.boundaryLeft);
	for (bit.First(); bit.Valid(); bit.Next()) {
		vort->Cell(bit) = - _v->dx_r(bit) + _u->dy_r(bit);
	}

    for(iit.First(); iit.Valid(); iit.Next()){ 
        // Calculating vorticity by du/dy - dv/dx
		vort->Cell(iit) = _u->dy_r(iit) - _v->dx_r(iit);
    }
    return vort;
}

/// Computes and returns the stream line values
const Grid * Compute::GetStream() {
    Grid * psi = new Grid(_geom, multi_real_t(1.0, 1.0));
    psi->Initialize(0.0);
    
    BoundaryIterator bit = BoundaryIterator(_geom);
    bit.SetBoundary(bit.boundaryLeft);
    bit.First();
    bit.Next();
    for (; bit.Valid(); bit.Next()) {
		psi->Cell(bit) = _geom->Mesh()[1]*_u->Cell(bit) + psi->Cell(bit.Down());
	}
    
    bit.SetBoundary(bit.boundaryBottom);
    bit.First();
    bit.Next();
    for (; bit.Valid(); bit.Next()) {
        psi->Cell(bit) = - _geom->Mesh()[0]*_v->Cell(bit) + psi->Cell(bit.Left());
	}
    
    InteriorIterator iit = InteriorIterator(_geom);
    
    for(iit.First(); iit.Valid(); iit.Next()){
        psi->Cell(iit) = - _geom->Mesh()[0]*_v->Cell(iit) + psi->Cell(iit.Left());
    }
    
    
	//psi->print();
	
    bit.SetBoundary(bit.boundaryLeft);
    bit.First();
    double offset = psi->Cell(bit);
    double buffer = 0.0;
    
    MPI_Status stat;
	
	for(int i = 0; i < _comm->ThreadDim()[1] - 1; i++) {
        if( _comm->getRank()==i*_comm->ThreadDim()[0]) {
            bit.SetBoundary(bit.boundaryTop);
            bit.First();
            buffer = offset + (double)psi->Cell(bit.Down());
	        MPI_Send(&buffer, 1, MPI_DOUBLE, (i+1)*_comm->ThreadDim()[0], 1, MPI_COMM_WORLD);
        }
        if( _comm->getRank()==(i+1)*_comm->ThreadDim()[0]) {
            MPI_Recv(&buffer, 1, MPI_DOUBLE, i*_comm->ThreadDim()[0], 1, MPI_COMM_WORLD, &stat);
            offset = buffer;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    for(int j = 0; j < _comm->ThreadDim()[0]-1; j++) {
        for(int i = 0; i < _comm->ThreadDim()[1]; i++) {
            if( _comm->getRank()==i*_comm->ThreadDim()[0]+j) {
                bit.SetBoundary(bit.boundaryRight);
                bit.First();
                buffer = offset + (double)psi->Cell(bit.Left());
                MPI_Send(&buffer, 1, MPI_DOUBLE, i*_comm->ThreadDim()[0]+j+1, i, MPI_COMM_WORLD);
            }
            if( _comm->getRank()==i*_comm->ThreadDim()[0]+j+1){
                MPI_Recv(&buffer, 1, MPI_DOUBLE, i*_comm->ThreadDim()[0]+j, i, MPI_COMM_WORLD, &stat);
                offset = buffer;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    psi->AddConstant(offset);

	return psi;
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
