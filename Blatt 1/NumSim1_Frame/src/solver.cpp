#include "solver.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include <math.h>

// Default Constructor
Solver::Solver() {

}

/// Constructor of the abstract Solver class
Solver::Solver(const Geometry * geom) {
    
    _geom = geom;
    
}

/// Destructor of the Solver Class
Solver::~Solver() {

}

/// Returns the residual at [it] for the pressure-Poisson equation
real_t Solver::localRes(const Iterator & it, const Grid * grid, const Grid * rhs) const {
    
    
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
    real_t h_x = _geom->Size()[0];
    real_t h_y = _geom->Size()[1];
    for(it.First(); it.Valid(); it.Next()){
        grid->Cell(it) = grid->Cell(it) + _omega*( rhs->Cell(it) -grid->dxx(it) - grid->dyy(it))/(-2.0/(h_x*h_x) - 2.0/(h_y*h_y));
    }
    
    
    real_t total_res = 0.0;
    for(it.First(); it.Valid(); it.Next()){
        total_res += localRes(it,grid,rhs)*localRes(it,grid,rhs);
    }
    total_res = total_res/(h_x*h_y);
    total_res = sqrt(total_res);
    return total_res;
    
    
}
