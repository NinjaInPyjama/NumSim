#include "geometry.hpp"

/// Constructs a default geometry:
// driven cavity with 128 x 128 grid, no-slip boundary conditions
// as shown below
//
//      u=1, v=0
//    -------------
//    |           |
// u=0|           |u=0
// v=0|           |v=0
//    |           |
//    |           |
//    -------------
//      u=0, v=0
Geometry::Geometry() {
    Load("default.geom");
	// Load("actual.geom");

	_h = multi_real_t(_length[0] / (_size[0] - 2), _length[1] / (_size[1] - 2));
}


Geometry::Geometry(const Communicator *comm) {
    Load("default.geom");
    _comm = comm;
    const multi_index_t real_size = multi_index_t((_bsize[0]-2)/_comm->ThreadDim()[0],(_bsize[1]-2)/_comm->ThreadDim()[1]) ;
    
    _redblack = (_comm->ThreadIdx()[0]*real_size[0]+_comm->ThreadIdx()[1]*real_size[1])%2 == 0;
    
    _size = multi_index_t(real_size[0]+2,real_size[1]+2);
    
    if ( _comm->isRight() ){
        _size[0] = _bsize[0] - (_comm->ThreadDim()[0]-1)*real_size[0];
    }
    if ( _comm->isTop() ){
        _size[1] = _bsize[1] - (_comm->ThreadDim()[1]-1)*real_size[1];
    }
    
    _h = multi_real_t(_blength[0] / (_bsize[0] - 2), _blength[1] / (_bsize[1] - 2));
    
    _length = multi_real_t(real_t(_size[0]-2)/(_bsize[0]-2)*_blength[0], real_t(_size[1]-2)/(_bsize[1]-2)*_blength[1]);
}



/// Loads a geometry from a file
void Geometry::Load(const char * file) {
    
    FILE* handle = fopen(file,"r");
    char name[20];
    multi_real_t inval_real;
    multi_index_t inval_index;
    while (!feof(handle)) {
        
        
        
        if (!fscanf(handle, "%s =", name)) continue;
        
        if (strcmp(name,"size") == 0) {
            if (fscanf(handle," %i %i\n",&inval_index[0],&inval_index[1])) {
                _bsize[0] = inval_index[0]+2;
                _bsize[1] = inval_index[1]+2;
            }
            continue;
        }
        
        if (strcmp(name,"length") == 0) {
            if (fscanf(handle," %lf %lf\n",&inval_real[0],&inval_real[1])) {
                _blength[0] = inval_real[0];
                _blength[1] = inval_real[1];
            }
            continue;
        }
        if (strcmp(name,"velocity") == 0) {
            if (fscanf(handle," %lf %lf\n",&inval_real[0],&inval_real[1])) {
                _velocity[0] = inval_real[0];
                _velocity[1] = inval_real[1];
            }
            continue;
        }
        if (strcmp(name,"pressure") == 0) {
            if (fscanf(handle," %lf\n",&inval_real[0]))
                _pressure = inval_real[0];
            continue;
        }

        if (strcmp(name,"geometry") == 0) {
            if (fscanf(handle," %s\n"),&name)
                if(strcmp(name,"free") == 0) {
					_field = new char[_bsize[1]][_bsize[0]];
                    for (int i = 0; i < _bsize[1]; i++) {
						fscanf(handle, " %s\n", _field[_bsize[1] - i - 1]);
                            //std::cout << line << std::endl;
                    }
                }
            continue;
            break;
        }
        
    }
	fclose(handle);
}

/// Returns whether the lower left corner is red or black
void Geometry::InitializeFlags(Grid * flag, Grid * type, Grid * value) const {
	value->Initialize(0.0);	

	Iterator it = Iterator(this);
    for(it.First(); it.Valid(); it.Next()) {
        multi_index_t pos = it.Pos();
        flag->Cell(it) = index_t(_field[pos[1]][pos[0]]);
    }

	index_t type_id = 0;
    for(it.First(); it.Valid(); it.Next()) {
		if (flag->Cell(it) != ' ') {
			type_id = 0;
			type_id += flag->Cell(it.Top()) == flag->Cell(it) ? 1 : 0;
			type_id += flag->Cell(it.Right()) == flag->Cell(it) ? 2 : 0;
			type_id += flag->Cell(it.Down()) == flag->Cell(it) ? 4 : 0;
			type_id += flag->Cell(it.Left()) == flag->Cell(it) ? 8 : 0;
			type->Cell(it) = type_id;
		}
		else type->Cell(it) = 0.0;
    }

    BoundaryIterator bit = BoundaryIterator(this);
    
    for(int i = 0 ; i < 4 ; i++) {
        bit.SetBoundary(i);
        int start_idx = -1;
        for(bit.First(); bit.Valid(); bit.Next()) {
            if(start_idx == -1 && (flag->Cell(it) == 'H' || flag->Cell(it) == 'V')) {
                start_idx = it.Value();
            }
            else if(start_idx != -1 && (flag->Cell(it) != 'H' || flag->Cell(it) != 'V')){
                int end_idx = it.Value();
				end_idx -= (i == bit.boundaryTop || i == bit.boundaryBottom) ? 1 : _size[0];

				//-4.0/(end_idx-start_idx+1.0)/(end_idx-start_idx+1.0)*velocity[1]*(y-start_idx+1.0/2.0)*(y-end_idx-1.0/2.0)
            
                BoundaryIterator bit_intern = BoundaryIterator(this);
                bit_intern.SetBoundary(i);
                for(bit_intern.First(); bit_intern.Valid(); bit_intern.Next()){
                    if(bit_intern.Value() <= end_idx && bit_intern.Value() >= start_idx){
                        value->Cell(bit) = -4.0/(end_idx-start_idx+1.0)/(end_idx-start_idx+1.0)*_velocity[1]*(bit.Value()-start_idx+1.0/2.0)*(bit.Value()-end_idx-1.0/2.0);
                    }
                }
                start_idx = -1;
            }
        }
    }
}


/// Returns whether the lower left corner is red or black
const bool & Geometry::RedBlack() const {
	return _redblack;
}

/// Returns the number of cells in each dimension
const multi_index_t & Geometry::Size() const {
	return _size;
}

/// Returns the total number of cells in each dimension
const multi_index_t & Geometry::TotalSize() const {
        return _bsize;
}

/// Returns the length of the domain
const multi_real_t & Geometry::Length() const {
	return _length;
}

/// Returns the total length of the domain
const multi_real_t & Geometry::TotalLength() const {
        return _blength;
}    

/// Returns the meshwidth
const multi_real_t & Geometry::Mesh() const {
	return _h;
}

/// Returns the initial velocity
const multi_real_t & Geometry::Velocity() const {
	return _velocity;
}

/// Returns the initial pressure
const real_t & Geometry::Pressure() const {
	return _pressure;
}


/// Updates the velocity field u
void Geometry::Update_U(Grid * u) const {
	// see script, p. 17
	BoundaryIterator bit = BoundaryIterator(this);
    
    // Iteration over right boundary
    if(_comm->isRight()) {
      bit.SetBoundary(bit.boundaryRight);
      for(bit.First(); bit.Valid(); bit.Next()) {
        u->Cell(bit) = 0.0;
        u->Cell(bit.Left()) = 0.0;
      }
    }
    
    // Iteration over left boundary
    if(_comm->isLeft()) {
      bit.SetBoundary(bit.boundaryLeft);
      for(bit.First(); bit.Valid(); bit.Next()) {
        u->Cell(bit) = 0.0;
      }
    }
    
    // Iteration over upper boundary
    if(_comm->isTop()) {
      bit.SetBoundary(bit.boundaryTop);
      for(bit.First(); bit.Valid(); bit.Next()) {
        u->Cell(bit) = 2.0 - u->Cell(bit.Down());
      }
    }
	
    // Iteration over lower boundary
    if(_comm->isBottom()){
      bit.SetBoundary(bit.boundaryBottom);
      for(bit.First(); bit.Valid(); bit.Next()) {
        u->Cell(bit) =  - u->Cell(bit.Top());
      }
    }
	
}

/// Updates the velocity field v
void Geometry::Update_V(Grid * v) const {
    // see script, p. 17	
    BoundaryIterator bit = BoundaryIterator(this);
    
    // Iteration over upper boundary
    if(_comm->isTop()) {
      bit.SetBoundary(bit.boundaryTop);
      for(bit.First(); bit.Valid(); bit.Next()) {
        v->Cell(bit) = 0.0;
        v->Cell(bit.Down()) = 0.0;
      }
    }
    
    // Iteration over lower boundary
    if(_comm->isBottom()) {
      bit.SetBoundary(bit.boundaryBottom);
      for(bit.First(); bit.Valid(); bit.Next()) {
        v->Cell(bit) = 0.0;
      }
    }
    
    // Iteration over left boundary
    if(_comm->isLeft()){
      bit.SetBoundary(bit.boundaryLeft);
      for(bit.First(); bit.Valid(); bit.Next()) {
        v->Cell(bit) = - v->Cell(bit.Right());
      }
    }
    
    // Iteration over right boundary
    if(_comm->isRight()){
      bit.SetBoundary(bit.boundaryRight);
      for(bit.First(); bit.Valid(); bit.Next()) {
        v->Cell(bit) = - v->Cell(bit.Left());
      }
    }
}

/// Updates the pressure field p
void Geometry::Update_P(Grid * p) const {
    // see script, p. 20 (p_{0,j} = p{1,j}), ...
    BoundaryIterator bit = BoundaryIterator(this);
    
    // Iteration over upper boundary
    if(_comm->isTop()){
      bit.SetBoundary(bit.boundaryTop);
      for(bit.First(); bit.Valid(); bit.Next()) {
			  p->Cell(bit) = p->Cell(bit.Down());
      }
    }
    
    // Iteration over right boundary
    if(_comm->isRight()){
      bit.SetBoundary(bit.boundaryRight);
      for(bit.First(); bit.Valid(); bit.Next()) {
        p->Cell(bit) = p->Cell(bit.Left());
      }
    }
    
    // Iteration over lower boundary
    if(_comm->isBottom()){
      bit.SetBoundary(bit.boundaryBottom);
      for(bit.First(); bit.Valid(); bit.Next()) {
        p->Cell(bit) = p->Cell(bit.Top());
      }
    }
    
    // Iteration over left boundary
    if(_comm->isLeft()){
      bit.SetBoundary(bit.boundaryLeft);
		  for(bit.First(); bit.Valid(); bit.Next()) {
			  p->Cell(bit) = p->Cell(bit.Right());
		  }
    }
}
