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
	return _tidx;
}

/** Returns the way the domain is partitioned among all processes
*/
const multi_index_t & Communicator::ThreadDim() const {
	return _tdim;
}

/** Returns whether this process is a red or a black field
*/
const bool & Communicator::EvenOdd() const {
	return _evenodd;
}

/** Gets the sum of all values and distributes the result among all
*  processes
*
* \param [in] val The data over which the sum is to be calculated
*/
real_t Communicator::gatherSum(const real_t & val) const {
	return real_t(0.0);
}

/** Finds the minimum of the values and distributes the result among
*  all processes
*
* \param [in] val The data over which to find the minimum
*/
real_t Communicator::gatherMin(const real_t & val) const {
	return real_t(0.0);
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
	return _tidx[0]==0;
}

/** Decide whether our right boundary is a domain boundary
*/
const bool Communicator::isRight() const {
	return _tidx[0]==_tdim[0]-1;
}

/** Decide whether our top boundary is a domain boundary
*/
const bool Communicator::isTop() const {
	return _tidx[1]==_tdim[1]-1;
}

/** Decide whether our bottom boundary is a domain boundary
*/
const bool Communicator::isBottom() const {
	return _tidx[1]==0;
}

/** Get MPI rank of current process
 */
const int & Communicator::getRank() const {
	return _rank;
}

/** Get number of MPI processes
*/
const int & Communicator::getSize() const {
	return _size;
}


  /** Function to sync ghost layer on left boundary:
   *  send values of own left boundary to left neighbor and
   *  and receive values from his right boundary
   *
   *   ------------ ------------
   *  |           x|y           |
   *  |           x|y           |
   *  |           x|y           |
   *  |           x|y           |
   *  |           x|y           |
   *   ------------ ------------
   *
   *   y: values that are sent
   *   x: values that are received
   *
   * \param [in] grid  values whose boundary shall be synced
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
