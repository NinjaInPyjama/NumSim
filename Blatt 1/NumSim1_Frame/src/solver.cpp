#include "solver.hpp"

// Default Constructor
Solver::Solver() {

}

/// Constructor of the abstract Solver class
Solver::Solver(const Geometry * geom) {

}

/// Destructor of the Solver Class
Solver::~Solver() {

}

/// Returns the residual at [it] for the pressure-Poisson equation
real_t Solver::localRes(const Iterator & it, const Grid * grid, const Grid * rhs) const {
	return 0;
}

/// Constructs an actual SOR solver
SOR::SOR(const Geometry * geom, const real_t & omega) {

}

/// Destructor
SOR::~SOR() {

}

/// Returns the total residual and executes a solver cycle
// @param grid current pressure values
// @param rhs right hand side
real_t SOR::Cycle(Grid * grid, const Grid * rhs) const {
	return 0;
}
