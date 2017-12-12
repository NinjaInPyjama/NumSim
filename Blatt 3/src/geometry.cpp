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
					_flag = new char[_bsize[1]*_bsize[0]];
                    for (int i = _bsize[1] - 1; i >= 0; i--) {
						char* line = new char[_bsize[0]];
						fscanf(handle, " %s\n", line);
						// std::cout << line << std::endl;
						for (int j = 0; j < _bsize[0]; j++) {
							_flag[i*_bsize[0] + j] = line[j];
						}
                    }
                }
            continue;
            break;
        }
        
    }
	fclose(handle);
}

void Geometry::InitializeFlags() const {

	Iterator it = Iterator(this);
	index_t type_id = 0;
    for(it.First(); it.Valid(); it.Next()) {
		if (_flag[it.Value()] != ' ') {
			type_id = 0;
			type_id += _flag[it.Top().Value()] == _flag[it.Value()] ? 1 : 0;
			type_id += _flag[it.Right().Value()] == _flag[it.Value()] ? 2 : 0;
			type_id += _flag[it.Down().Value()] == _flag[it.Value()] ? 4 : 0;
			type_id += _flag[it.Left().Value()] == _flag[it.Value()] ? 8 : 0;
			_type[it.Value()] = type_id;
		}
		else _type[it.Value()] = 0.0;
    }

    BoundaryIterator bit = BoundaryIterator(this);
    
    for(int i = 0 ; i < 4 ; i++) {
        bit.SetBoundary(i);
        int start_idx = -1; // saves first position of 'V' or 'H'
        for(bit.First(); bit.Valid(); bit.Next()) {
            if(start_idx == -1 && (_flag[bit.Value()] == 'H' || _flag[bit.Value()] == 'V')) {
                start_idx = bit.Value();
            }
            else if(start_idx != -1 && (_flag[bit.Value()] != 'H' || _flag[bit.Value()] != 'V')){
                int end_idx = bit.Value();
				end_idx -= (i == bit.boundaryTop || i == bit.boundaryBottom) ? 1 : _size[0];

				//-4.0/(end_idx-start_idx+1.0)/(end_idx-start_idx+1.0)*velocity[1]*(y-start_idx+1.0/2.0)*(y-end_idx-1.0/2.0)
            
                BoundaryIterator bit_intern = BoundaryIterator(this);
                bit_intern.SetBoundary(i);
                for(bit_intern.First(); bit_intern.Valid(); bit_intern.Next()){
                    if(bit_intern.Value() <= end_idx && bit_intern.Value() >= start_idx){
                        _value[bit_intern.Value()] = -4.0/(end_idx-start_idx+1.0)/(end_idx-start_idx+1.0)*_velocity[(i+1)%2]*(bit.Value()-start_idx+1.0/2.0)*(bit.Value()-end_idx-1.0/2.0);
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

	Iterator it(this);
	for (it.First(); it.Valid(); it.Next()) {
		switch (_flag[it.Value()]) {
		case '#': // NOSLIP
			switch (_type[it.Value()]) {
			case boundaryTop:
				v->Cell(it) = 0.0;
				v->Cell(it.Down()) = 0.0;
				break;
			case boundaryRight:
				v->Cell(it) = -v->Cell(it.Left());
				break;
			case boundaryBottom:
				v->Cell(it) = 0.0;
				break;
			case boundaryLeft:
				v->Cell(it) = -v->Cell(it.Right());
				break;
			case cornerTopRight:
				v->Cell(it.Down()) = 0.0;
				v->Cell(it) = -v->Cell(it.Left());
				break;
			case cornerBottomRight:
				v->Cell(it) = 0.0;
				break;
			case cornerBottomLeft:
				v->Cell(it) = 0.0;
				break;
			case cornerTopLeft:
				v->Cell(it.Down()) = 0.0;
				v->Cell(it) = -v->Cell(it.Right());
				break;
			case inner:
				v->Cell(it) = 0.0;
				break;
			default:
				break;
			}
			break;
		case '-': // Horizontal SLIP
			switch (_type[it.Value()]) {
			case boundaryTop:
				v->Cell(it) = 0.0;
				v->Cell(it.Down()) = 0.0;
				break;
			case boundaryRight:
				v->Cell(it) = -v->Cell(it.Left());
				break;
			case boundaryBottom:
				v->Cell(it) = 0.0;
				break;
			case boundaryLeft:
				v->Cell(it) = -v->Cell(it.Right());
				break;
			case cornerTopRight:
				v->Cell(it.Down()) = 0.0;
				v->Cell(it) = -v->Cell(it.Left());
				break;
			case cornerBottomRight:
				v->Cell(it) = 0.0;
				break;
			case cornerBottomLeft:
				v->Cell(it) = 0.0;
				break;
			case cornerTopLeft:
				v->Cell(it.Down()) = 0.0;
				v->Cell(it) = -v->Cell(it.Right());
				break;
			case inner:
				v->Cell(it) = 0.0;
				break;
			default:
				break;
			}
			break;
		case '|': // Vertical SLIP
			switch (_type[it.Value()]) {
			case boundaryTop:
				v->Cell(it) = v->Cell(it.Left());
				v->Cell(it.Down()) = v->Cell(it.Down().Left());
				break;
			case boundaryRight:
				v->Cell(it) = v->Cell(it.Left());
				break;
			case boundaryBottom:
				v->Cell(it) = v->Cell(it.Left());
				break;
			case boundaryLeft:
				v->Cell(it) = v->Cell(it.Right());
				break;
			case cornerTopRight:
				v->Cell(it) = v->Cell(it.Left());
				v->Cell(it.Down()) = v->Cell(it.Down().Left());
				break;
			case cornerBottomRight:
				v->Cell(it) = v->Cell(it.Left());
				break;
			case cornerBottomLeft:
				v->Cell(it) = v->Cell(it.Right());
				break;
			case cornerTopLeft:
				v->Cell(it) = v->Cell(it.Right());
				v->Cell(it.Down()) = v->Cell(it.Down().Right());
				break;
			case inner:
				v->Cell(it) = 0.0;
				break;
			default:
				break;
			}
		case 'O': // OUTFLOW
			switch (_type[it.Value()]) {
			case boundaryTop:
				v->Cell(it.Down()) = v->Cell(it.Down().Down());
				v->Cell(it) = v->Cell(it.Down());
				break;
			case boundaryRight:
				v->Cell(it) = v->Cell(it.Left());
				break;
			case boundaryBottom:
				v->Cell(it) = v->Cell(it.Top());
				break;
			case boundaryLeft:
				v->Cell(it) = v->Cell(it.Right());
				break;
			case cornerTopRight:
				break;
			case cornerBottomRight:
				break;
			case cornerBottomLeft:
				break;
			case cornerTopLeft:
				break;
			case inner:
				break;
			default:
				break;
			}
		case 'V': // Vertical INFLOW
			break;
		case 'H': // Horizontal INFLOW
			v->Cell(it) = _value[it.Value()];
			break;
		default:
			break;
		}
	}
}

/// Updates the pressure field p
void Geometry::Update_P(Grid * p) const {
    // see script, p. 20 (p_{0,j} = p{1,j}), ...
	Iterator it(this);
	for (it.First(); it.Valid(); it.Next()) {
		switch (_flag[it.Value()]) {
		case '#': // NOSLIP
			switch (_type[it.Value()]) {
			case boundaryTop:
				p->Cell(it) = p->Cell(it.Down());
				break;
			case boundaryRight:
				p->Cell(it) = p->Cell(it.Left());
				break;
			case boundaryBottom:
				p->Cell(it) = p->Cell(it.Top());
				break;
			case boundaryLeft:
				p->Cell(it) = p->Cell(it.Right());
				break;
			case cornerTopRight:
				p->Cell(it) = 0.5*(p->Cell(it.Left()) + p->Cell(it.Down()));
				break;
			case cornerBottomRight:
				p->Cell(it) = 0.5*(p->Cell(it.Left()) + p->Cell(it.Top()));
				break;
			case cornerBottomLeft:
				p->Cell(it) = 0.5*(p->Cell(it.Right()) + p->Cell(it.Top()));
				break;
			case cornerTopLeft:
				p->Cell(it) = 0.5*(p->Cell(it.Right()) + p->Cell(it.Down()));
				break;
			case inner:
				p->Cell(it) = 0.0;
				break;
			default:
				break;
			}
			break;
		case '-': // Horizontal SLIP
			switch (_type[it.Value()]) {
			case boundaryTop:
				p->Cell(it) = p->Cell(it.Down());
				break;
			case boundaryRight:
				p->Cell(it) = 2 * _pressure - p->Cell(it.Left());
				break;
			case boundaryBottom:
				p->Cell(it) = p->Cell(it.Top());
				break;
			case boundaryLeft:
				p->Cell(it) = 2 * _pressure - p->Cell(it.Right());
				break;
			case cornerTopRight:
				p->Cell(it) = 0.5*(p->Cell(it.Left()) + p->Cell(it.Down()));
				break;
			case cornerBottomRight:
				p->Cell(it) = 0.5*(p->Cell(it.Left()) + p->Cell(it.Top()));
				break;
			case cornerBottomLeft:
				p->Cell(it) = 0.5*(p->Cell(it.Right()) + p->Cell(it.Top()));
				break;
			case cornerTopLeft:
				p->Cell(it) = 0.5*(p->Cell(it.Right()) + p->Cell(it.Down()));
				break;
			case inner:
				p->Cell(it) = 0.0;
				break;
			default:
				break;
			}
			break;
		case '|': // Vertical SLIP
			switch (_type[it.Value()]) {
			case boundaryTop:
				p->Cell(it) = 2 * _pressure - p->Cell(it.Down());
				break;
			case boundaryRight:
				p->Cell(it) = p->Cell(it.Left());
				break;
			case boundaryBottom:
				p->Cell(it) = 2 * _pressure - p->Cell(it.Top());
				break;
			case boundaryLeft:
				p->Cell(it) = p->Cell(it.Right());
				break;
			case cornerTopRight:
				p->Cell(it) = 0.5*(p->Cell(it.Left()) + p->Cell(it.Down()));
				break;
			case cornerBottomRight:
				p->Cell(it) = 0.5*(p->Cell(it.Left()) + p->Cell(it.Top()));
				break;
			case cornerBottomLeft:
				p->Cell(it) = 0.5*(p->Cell(it.Right()) + p->Cell(it.Top()));
				break;
			case cornerTopLeft:
				p->Cell(it) = 0.5*(p->Cell(it.Right()) + p->Cell(it.Down()));
				break;
			case inner:
				p->Cell(it) = 0.0;
				break;
			default:
				break;
			}
		case 'O': // OUTFLOW
			p->Cell(it) = 0.0;
		case 'V': // Vertical INFLOW
			break;
		case 'H': // Horizontal INFLOW
			break;
		default:
			break;
		}
	}
}
