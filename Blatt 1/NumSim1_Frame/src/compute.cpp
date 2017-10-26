#include "compute.hpp"

/// Creates a compute instance with given geometry and parameter
Compute::Compute(const Geometry * geom, const Parameter * param) {

}

/// Deletes all grids
Compute::~Compute() {

}


/// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual
// etc.)
void Compute::TimeStep(bool printInfo) {

}

/// Returns the simulated time in total
const real_t & Compute::GetTime() const {
	return 0;
}


/// Returns the pointer to U
const Grid * Compute::GetU() const {
	return nullptr;
}

/// Returns the pointer to V
const Grid * Compute::GetV() const {
	return nullptr;
}

/// Returns the pointer to P
const Grid * Compute::GetP() const {
	return nullptr;
}

/// Returns the pointer to RHS
const Grid * Compute::GetRHS() const {
	return nullptr;
}


/// Computes and returns the absolute velocity
const Grid * Compute::GetVelocity() {
	return nullptr;
}

/// Computes and returns the vorticity
const Grid * Compute::GetVorticity() {
	return nullptr;
}

/// Computes and returns the stream line values
const Grid * Compute::GetStream() {
	return nullptr;
}


/// Compute the new velocites u,v
void Compute::NewVelocities(const real_t & dt) {

}

/// Compute the temporary velocites F,G
void Compute::MomentumEqu(const real_t & dt) {

}

/// Compute the RHS of the poisson equation
void Compute::RHS(const real_t & dt) {

}