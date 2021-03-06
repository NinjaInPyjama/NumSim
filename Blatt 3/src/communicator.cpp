#include "communicator.hpp"
#include <mpi/mpi.h>

/** Communicator constructor; initializes MPI Environment
*
* \param [in] argc Number of arguments program was started with
* \param [in] argv Arguments passed to the program on start
*/
Communicator::Communicator(int * argc, char *** argv) {
    MPI_Init(argc, argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    
    MPI_Comm_size(MPI_COMM_WORLD, &this->_size);   
    
    for (int i=(int)sqrt(_size) ; i>=1;i--){
        if((real_t( _size)/real_t(i)) == _size/i){
            _tdim = multi_index_t(i,_size/i);
            break;
        }
    }
    //unten links ist (0,0)
    _tidx[0] = _rank%_tdim[0];
    _tidx[1] = (index_t)(_rank/_tdim[0]);
    
    //unten links ist even (true)
    _evenodd = (_tidx[0]+_tidx[1])%2 == 0;
}


/** Communicator destructor; finalizes MPI Environment
*/
Communicator::~Communicator() {
    MPI_Finalize();
}


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
	real_t sum = 0;
    MPI_Allreduce(&val, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return sum;
}

/** Finds the minimum of the values and distributes the result among
*  all processes
*
* \param [in] val The data over which to find the minimum
*/
real_t Communicator::gatherMin(const real_t & val) const {
	real_t min = 0;
    MPI_Allreduce(&val, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    
    return min;
}

/** Finds the maximum of the values and distributes the result among
*  all processes
*
* \param [in] val The data over which to find the maximum
*/
real_t Communicator::gatherMax(const real_t & val) const {
	real_t max = 0;
    MPI_Allreduce(&val, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    return max;
}

/** Synchronizes ghost layer
*
* \param [in] grid  The values to sync
*/
void Communicator::copyBoundary(Grid * grid) const {
	
    if(this->EvenOdd()) {
        if(!this->isLeft()) {
            this->copyLeftBoundary(grid);
        }
        if(!this->isRight()) {
            this->copyRightBoundary(grid);
        }
        if(!this->isTop()) {
            this->copyTopBoundary(grid);
        }
        if(!this->isBottom()) {
            this->copyBottomBoundary(grid);
        }
    } else {
        if(!this->isRight()) {
            this->copyRightBoundary(grid);
        }
        if(!this->isLeft()) {
            this->copyLeftBoundary(grid);
        }
        if(!this->isBottom()) {
            this->copyBottomBoundary(grid);
        }
        if(!this->isTop()) {
            this->copyTopBoundary(grid);
        }
    }
}

/** Decide whether our left boundary is a domain boundary
*/
const bool Communicator::isLeft() const {
	return _tidx[0]==0;
}

/** Decide whether our right boundary is a domain boundary
*/
const bool Communicator::isRight() const {
	return _tidx[0]==(_tdim[0]-1);
}

/** Decide whether our top boundary is a domain boundary
*/
const bool Communicator::isTop() const {
	return _tidx[1]==(_tdim[1]-1);
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
	MPI_Status stat;
	const index_t bufferlength = grid->getGeometry()->Size()[1];
	real_t buffer[bufferlength];
	BoundaryIterator bit = BoundaryIterator(grid->getGeometry());
	bit.SetBoundary(bit.boundaryLeft);

    int i = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		buffer[i] = grid->Cell(bit.Right());
		i++;
	}
	
	int result = MPI_Sendrecv_replace(buffer, bufferlength, MPI_DOUBLE, _rank - 1, 0, _rank - 1, 1, MPI_COMM_WORLD, &stat);

	i = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		grid->Cell(bit) = buffer[i];
		i++;
	}
    
    return (result == MPI_SUCCESS);
}

/** Function to sync ghost layer on right boundary
*  Details analog to left boundary
*
* \param [in] grid  values whose boundary shall be synced
*/
bool Communicator::copyRightBoundary(Grid * grid) const {
	MPI_Status stat;
	const index_t bufferlength = grid->getGeometry()->Size()[1];
	real_t buffer[bufferlength];
	BoundaryIterator bit = BoundaryIterator(grid->getGeometry());
	bit.SetBoundary(bit.boundaryRight);

    int i = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		buffer[i] = grid->Cell(bit.Left());
		i++;
	}    
	
	int result = MPI_Sendrecv_replace(buffer, bufferlength, MPI_DOUBLE, _rank + 1, 1, _rank + 1, 0, MPI_COMM_WORLD, &stat);

	i = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		grid->Cell(bit) = buffer[i];
		i++;
	}
	
	return (result == MPI_SUCCESS);
}

/** Function to sync ghost layer on top boundary
*  Details analog to left boundary
*
* \param [in] grid  values whose boundary shall be synced
*/
bool Communicator::copyTopBoundary(Grid * grid) const {
	MPI_Status stat;
	const index_t bufferlength = grid->getGeometry()->Size()[0];
	real_t buffer[bufferlength];
	BoundaryIterator bit = BoundaryIterator(grid->getGeometry());
	bit.SetBoundary(bit.boundaryTop);

    int i = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		buffer[i] = grid->Cell(bit.Down());
		i++;
	}
    
	int result = MPI_Sendrecv_replace(buffer, bufferlength, MPI_DOUBLE, _rank + _tdim[0], 3, _rank + _tdim[0], 2, MPI_COMM_WORLD, &stat);
    
	i = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		grid->Cell(bit) = buffer[i];
		i++;
	}
	
	return (result == MPI_SUCCESS);
}

/** Function to sync ghost layer on bottom boundary
*  Details analog to left boundary
*
* \param [in] grid  values whose boundary shall be synced
/** Get
*/
bool Communicator::copyBottomBoundary(Grid * grid) const {
	MPI_Status stat;
	const index_t bufferlength = grid->getGeometry()->Size()[0];
	real_t buffer[bufferlength];
	BoundaryIterator bit = BoundaryIterator(grid->getGeometry());
    bit.SetBoundary(bit.boundaryBottom);

    int i = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		buffer[i] = grid->Cell(bit.Top());
		i++;
	}
    
	int result = MPI_Sendrecv_replace(buffer, bufferlength, MPI_DOUBLE, _rank - _tdim[0], 2, _rank - _tdim[0], 3, MPI_COMM_WORLD, &stat);

    i = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		grid->Cell(bit) = buffer[i];
		i++;
	}
	
	return (result == MPI_SUCCESS);
}
