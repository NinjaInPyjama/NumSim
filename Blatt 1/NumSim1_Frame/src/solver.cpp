#include "solver.hpp"

// Default Constructor
Solver::Solver() {}

/// Constructor of the abstract Solver class
Solver::Solver(const Geometry * geom) {
    _geom = geom;
}

/// Destructor of the Solver Class
Solver::~Solver() {

}

/// Returns the residual at [it] for the pressure-Poisson equation
real_t Solver::localRes(const Iterator & it, const Grid * grid, const Grid * rhs) const {
    // see script, p. 25
	// TODO: alternativ? => 
	// return grid->Cell(it) - (dx * dx * dy * dy) / (2.0 * (dx * dx + dy * dy)) * ((grid->Cell(it.Left()) + grid->Cell(it.Right()))/(dx*dx) + (grid->Cell(it.Down()) + grid->Cell(it.Top())) / (dy*dy) - rhs->Cell(it));
	return grid->dxx(it) + grid->dyy(it) - rhs->Cell(it);
}

/// Constructs an actual SOR solver
SOR::SOR(const Geometry * geom, const real_t & omega) {
    _omega = omega;
    _geom = geom;
}

/// Constructs an actual SOR solver
SOR::SOR(const Geometry * geom) {
    _omega = real_t(2.0/(1.0 + sin(M_PI*geom->Mesh()[0])));
    _geom = geom;
}

/// Destructor
SOR::~SOR() {
    
}

/// Returns the total residual and executes a solver cycle
// @param grid current pressure values
// @param rhs right hand side
real_t SOR::Cycle(Grid * grid, const Grid * rhs) const {
    InteriorIterator it = InteriorIterator(_geom);
    real_t dx = _geom->Mesh()[0];
    real_t dy = _geom->Mesh()[1];

    for(it.First(); it.Valid(); it.Next()){
		// TODO: alternativ?! =>
		// grid->Cell(it) = (1 - _omega) * grid->Cell(it) + _omega * (dx * dx * dy * dy) / (2.0 * (dx * dx + dy * dy)) * ((grid->Cell(it.Left()) + grid->Cell(it.Right()))/(dx*dx) + (grid->Cell(it.Down()) + grid->Cell(it.Top())) / (dy*dy) - rhs->Cell(it));
		// see script, p. 26, formular (4.1)
        grid->Cell(it) = grid->Cell(it) + _omega*( rhs->Cell(it) -grid->dxx(it) - grid->dyy(it))/(-2.0/(dx*dx) - 2.0/(dy*dy));
    }
    
	_geom->Update_P(grid);

    real_t total_res = 0.0;
	real_t local_res = 0.0;

    for(it.First(); it.Valid(); it.Next()){
		local_res = localRes(it, grid, rhs);
		total_res += local_res*local_res;
    }

    total_res = total_res/(_geom->Size()[0] * _geom->Size()[1]);
    total_res = sqrt(total_res);
    return total_res;
}
