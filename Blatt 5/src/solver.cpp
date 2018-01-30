#include "solver.hpp"

// Default Constructor
Solver::Solver() {}

/// Constructor of the abstract Solver class
Solver::Solver(const Geometry * geom) {
    _geom = geom;
}

/// Destructor of the Solver Class
Solver::~Solver() {}

/// Returns the residual at [it] for the pressure-Poisson equation
real_t Solver::localRes(const Iterator & it, const Grid * grid, const Grid * rhs) const {
    // see script, p. 25
	return grid->dxx(it) + grid->dyy(it) - rhs->Cell(it);
}

/// Constructs an actual SOR solver
SOR::SOR(const Geometry * geom, const real_t & omega) {
    _omega = omega;
    _geom = geom;
}

/// Constructs an actual SOR solver
SOR::SOR(const Geometry * geom) {
    _omega = std::min(real_t(2.0/(1.0 + sin(M_PI*geom->Mesh()[0]))), real_t(2.0/(1.0 + sin(M_PI*geom->Mesh()[1])))); // optimal omega, see script p. 26
    _geom = geom;
}

/// Destructor
SOR::~SOR() {}

/// Returns the total residual and executes a solver cycle
// @param grid current pressure values
// @param rhs right hand side
real_t SOR::Cycle(Grid * grid, const Grid * rhs) const {
    InteriorIterator it = InteriorIterator(_geom);
    real_t dx = _geom->Mesh()[0];
    real_t dy = _geom->Mesh()[1];
    real_t factor = 1.0/(2.0/(dx*dx) + 2.0/(dy*dy));
    real_t total_res = 0.0;
	real_t local_res = 0.0;
    index_t fluidcells = 0;

    for(it.First(); it.Valid(); it.Next()){
		// see script, p. 26, formular (4.1)
        
        fluidcells++;
        local_res = localRes(it, grid, rhs);
        total_res += local_res*local_res;
        grid->Cell(it) = grid->Cell(it) + _omega*local_res*factor;
    }
    
	//Updating boundary values to reduce artifactial values in calculation of residual
	_geom->Update_P(grid);

    total_res = total_res/real_t(fluidcells);
    return total_res;
}


MGSolver::MGSolver(const Geometry *geom, index_t smoothCycles) {
    _geom = geom;
    _smoothCycles = smoothCycles;
}


MGSolver::MGSolver(const Geometry *geom) {
    _geom = geom;
    _smoothCycles = 3;
}
  
/// Destructor
MGSolver::~MGSolver() {}

void MGSolver::MGCycle(Grid *grid, Grid *rhs) {
    MGCycle(new MultiGrid(grid, _geom), new MultiGrid(rhs, _geom), 1);
}
  
real_t MGSolver::localRes(const MGIterator &it, MultiGrid *grid, MultiGrid *rhs) {
    return grid->dxx(it) + grid->dyy(it) - rhs->Cell(it);
}
  
void MGSolver::smooth(MultiGrid *grid, MultiGrid *rhs) {
    
}
  
MultiGrid * MGSolver::RES( MultiGrid *grid, MultiGrid *rhs) {
    MultiGrid * res = new MultiGrid(grid, _geom);
    MGInteriorIterator it = MGInteriorIterator(_geom, res);
    for(it.First(); it.Valid(); it.Next()) {
      res->Cell(it) = -rhs->Cell(it);
    }
    return res;
}

void MGSolver::MGCycle(MultiGrid *grid, MultiGrid *rhs, index_t cellSize) {
    for(int i=0; i<_smoothCycles; i++) smooth(grid, rhs); // pre-smoothing
    if(cellSize < 0.5*std::min(_geom->Size()[0], _geom->Size()[1])) {
      // TODO: solver Cycle
    }
    for(int i=0; i<_smoothCycles; i++) smooth(grid, rhs); // post-smoothing
}
