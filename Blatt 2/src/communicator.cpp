#include "communicator.hpp"

/** Communicator constructor; initializes MPI Environment
*
* \param [in] argc Number of arguments program was started with
* \param [in] argv Arguments passed to the program on start
*/
Communicator::Communicator(int * argc, char *** argv) {

}


/** Communicator destructor; finalizes MPI Environment
*/
Communicator::~Communicator() {}


/** Returns the position of the current process with respect to the
*  fields lower left corner
*/
const multi_index_t & Communicator::ThreadIdx() const {
	return multi_index_t();
}

/** Returns the way the domain is partitioned among all processes
*/
const multi_index_t & Communicator::ThreadDim() const {
	return multi_index_t();
}

/** Returns whether this process is a red or a black field
*/
const bool & Communicator::EvenOdd() const {
	return false;
}

/** Gets the sum of all values and distributes the result among all
*  processes
*
* \param [in] val The data over which the sum is to be calculated
*/
real_t Communicator::gatherSum(const real_t & val) const {
	return real_t();
}

/** Finds the minimum of the values and distributes the result among
*  all processes
*
* \param [in] val The data over which to find the minimum
*/
real_t Communicator::gatherMin(const real_t & val) const {
	return real_t();
}

/** Finds the maximum of the values and distributes the result among
*  all processes
*
* \param [in] val The data over which to find the maximum
*/
real_t Communicator::gatherMax(const real_t & val) const {
	return real_t();
}

/** Synchronizes ghost layer
*
* \param [in] grid  The values to sync
*/
void Communicator::copyBoundary(Grid * grid) const {

}

/** Decide whether our left boundary is a domain boundary
*/
const bool Communicator::isLeft() const {
	return false;
}

/** Decide whether our right boundary is a domain boundary
*/
const bool Communicator::isRight() const {
	return false;
}

/** Decide whether our top boundary is a domain boundary
*/
const bool Communicator::isTop() const {
	return false;
}

/** Decide whether our bottom boundary is a domain boundary
*/
const bool Communicator::isBottom() const {
	return false;
}

/** Decide whether our bottom boundary is a domain boundary
*/
const int & Communicator::getRank() const {
	// TODO: insert return statement here
}

/** Get number of MPI processes
*/
const int & Communicator::getSize() const {
	// TODO: insert return statement here
}


/** Get number of MPI processes
*/
bool Communicator::copyLeftBoundary(Grid * grid) const {
	return false;
}

/** Function to sync ghost layer on right boundary
*  Details analog to left boundary
*
* \param [in] grid  values whose boundary shall be synced
*/
bool Communicator::copyRightBoundary(Grid * grid) const {
	return false;
}

/** Function to sync ghost layer on top boundary
*  Details analog to left boundary
*
* \param [in] grid  values whose boundary shall be synced
*/
bool Communicator::copyTopBoundary(Grid * grid) const {
	return false;
}

/** Function to sync ghost layer on bottom boundary
*  Details analog to left boundary
*
* \param [in] grid  values whose boundary shall be synced
*/
bool Communicator::copyBottomBoundary(Grid * grid) const {
	return false;
}
