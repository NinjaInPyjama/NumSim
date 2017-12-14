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
    
    std::cout << "Hallo3.1" << std::endl;
    Load("default.geom");
    _comm = comm;
    std::cout << "Halloooo" << std::endl;
    const multi_index_t real_size = multi_index_t((_bsize[0])/_comm->ThreadDim()[0],(_bsize[1])/_comm->ThreadDim()[1]) ;
    
    _redblack = (_comm->ThreadIdx()[0]*real_size[0]+_comm->ThreadIdx()[1]*real_size[1])%2 == 0;
    
    _size = multi_index_t(real_size[0],real_size[1]);
    
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
    std::cout << "Hallo3.2" << std::endl;
    while (!feof(handle)) {
        
        
        
        if (!fscanf(handle, "%s =", name)) continue;
        
        if (strcmp(name,"size") == 0) {
            if (fscanf(handle," %i %i\n",&inval_index[0],&inval_index[1])) {
                _bsize[0] = inval_index[0];
                _bsize[1] = inval_index[1];
            }
            std::cout << "Hallo3.3" << std::endl;
            continue;
        }
        
        if (strcmp(name,"length") == 0) {
            if (fscanf(handle," %lf %lf\n",&inval_real[0],&inval_real[1])) {
                _blength[0] = inval_real[0];
                _blength[1] = inval_real[1];
            }
            std::cout << "Hallo3.4" << std::endl;
            continue;
        }
        if (strcmp(name,"velocity") == 0) {
            if (fscanf(handle," %lf %lf\n",&inval_real[0],&inval_real[1])) {
                _velocity[0] = inval_real[0];
                _velocity[1] = inval_real[1];
                std::cout << "Hallo3.5" << std::endl;
            }
            continue;
        }
        if (strcmp(name,"pressure") == 0) {
            if (fscanf(handle," %lf\n",&inval_real[0]))
                _pressure = inval_real[0];
            std::cout << "Hallo3.6" << std::endl;
            continue;
        }

        if (strcmp(name,"geometry") == 0) {
            if (fscanf(handle," %s\n", name))
                std::cout << "Hallo3.7  " << name << std::endl;
                if(strcmp(name,"free") == 0) {
                    
                    std::cout << "Hallo3.8" << std::endl;
                    
					_flag = new char[_bsize[1]*_bsize[0]];
                    _type = new int[_bsize[1]*_bsize[0]];
                    _value = new real_t[_bsize[1]*_bsize[0]];
                    
                    char* line = new char[_bsize[0]];
                    std::cout << "Hallo3.8.1 " << _flag[1] << "k" << std::endl;
                    for (int i = _bsize[1] - 1; i >= 0; i--) {
						
//                         for(int k = 0; k<_bsize[0]-1; k++) {
//                             fscanf(handle, "%c[1]", &line[k]);
//                             //std::cout << line[k] << std::endl;
//                         }
//                         fscanf(handle, "%c[1]\n", &line[_bsize[0]-1]);
                        fscanf(handle, "%[-HV| #IO]\n", line);
                        std::cout << "Hallo3.8.2 " << i << line << std::endl;
                        //std::cout << line << std::endl;
						for (int j = 0; j < _bsize[0]; j++) {
							_flag[i*_bsize[0] + j] = line[j];
                            _type [i*_bsize[0] + j] = 0;
                            _value [i*_bsize[0] + j] = 0.0;
						}
						
                    }
                    std::cout << "Hallo3.8.3 " << _flag[160-1] << std::endl;
                }
            continue;
            break;
        }
        
    }
	fclose(handle);
}

void Geometry::InitializeFlags() const {

	Iterator it = Iterator(this);
    std::cout << "Hallo4.1" << std::endl;
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
        int end_idx = -1;
        for(bit.First(); bit.Valid(); bit.Next()) {
            if(start_idx == -1 && (_flag[bit.Value()] == 'H' || _flag[bit.Value()] == 'V')) {
                start_idx = bit.Value();
            }
            else if(start_idx != -1 && (_flag[bit.Value()] != 'H' && _flag[bit.Value()] != 'V')){
                end_idx = bit.Value();
				end_idx -= (i == bit.boundaryTop || i == bit.boundaryBottom) ? 1 : -_size[0];

				//-4.0/(end_idx-start_idx+1.0)/(end_idx-start_idx+1.0)*velocity[1]*(y-start_idx+1.0/2.0)*(y-end_idx-1.0/2.0)
            
                BoundaryIterator bit_intern = BoundaryIterator(this);
                bit_intern.SetBoundary(i);
                std::cout << start_idx << "  " << end_idx << std::endl;
                
                for(bit_intern.First(); bit_intern.Valid(); bit_intern.Next()){
                    if(bit_intern.Value() <= end_idx && bit_intern.Value() >= start_idx){
                        _value[bit_intern.Value()] = -4.0/((end_idx-start_idx+1.0)*(end_idx-start_idx+1.0))*_velocity[(i+1)%2]*(real_t(bit_intern.Value())-real_t(start_idx)+1.0/2.0)*(real_t(bit_intern.Value())-real_t(end_idx)-1.0/2.0);
                    }
                }
                start_idx = -1;
                end_idx = -1;
            }
        }
    }
    
}


/// Prints the values of grid
void Geometry::print(const char & c) const {
	std::cout.precision(2);
    if(c == 'v') {
	for (int i = _size[1] - 1; i >= 0; i--) {
		for (int j = 0; j < _size[0]; j++) {
			std::cout  << " " << _value[i*_size[0] + j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
    }else if(c == 'f') {
    for (int i = _size[1] - 1; i >= 0; i--) {
		for (int j = 0; j < _size[0]; j++) {
			std::cout  << _flag[i*_size[0] + j] ;
		}
		std::cout << std::endl;
	}
    }else if(c == 't') {
    for (int i = _size[1] - 1; i >= 0; i--) {
		for (int j = 0; j < _size[0]; j++) {
			std::cout <<  _type[i*_size[0] + j] ;
		}
		std::cout << std::endl;
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
	Iterator it = Iterator(this);

	for (it.First(); it.Valid(); it.Next()) {
		switch (_flag[it.Value()]) {
		case '#': // NOSLIP
			switch (_type[it.Value()]) {
			case boundaryBottom:
				u->Cell(it) = - u->Cell(it.Top());
				break;
			case boundaryTop:
				u->Cell(it) = - u->Cell(it.Down());
				break;
			case boundaryRight:
				u->Cell(it) = 0.0;
				u->Cell(it.Left()) = 0.0;
				break;
			case boundaryLeft:
				u->Cell(it) = 0.0;
				break;
			case cornerTopRight:
				u->Cell(it) = - u->Cell(it.Down());
				u->Cell(it.Left()) = 0.0;
				break;
			case cornerTopLeft:
				u->Cell(it) = 0.0;
				break;
			case cornerBottomLeft:
				u->Cell(it) = 0.0;
				break;
			case cornerBottomRight:
				u->Cell(it) = - u->Cell(it.Top());
				u->Cell(it.Left()) = 0;
				break;
			case inner:
				u->Cell(it) = 0.0;
				break;
			default:
				break;
			}
			break;
		case '-': // Horizontal SLIP 
			switch (_type[it.Value()]) {
			case boundaryBottom:
				u->Cell(it) = u->Cell(it.Top());
				break;
			case boundaryTop:
				u->Cell(it) = u->Cell(it.Down());
				break;
			case boundaryRight:
				u->Cell(it.Left()) = u->Cell(it.Left().Left());
				u->Cell(it) = u->Cell(it.Left());
				break;
			case boundaryLeft:
				u->Cell(it) = u->Cell(it.Right());
				break;
			case cornerTopRight:
				u->Cell(it.Left()) = u->Cell(it.Left().Left());
				u->Cell(it) = u->Cell(it.Down());
				break;
			case cornerTopLeft:
				u->Cell(it) = u->Cell(it.Right());
				break;
			case cornerBottomLeft:
				u->Cell(it) = u->Cell(it.Right());
				break;
			case cornerBottomRight:
				u->Cell(it.Left()) = u->Cell(it.Left().Left());
				u->Cell(it) = u->Cell(it.Top());
				break;
			}
			break;
		case '|': // Vertical SLIP 
			switch (_type[it.Value()]) {
			case boundaryBottom:
				u->Cell(it) = - u->Cell(it.Top());
				break;
			case boundaryTop:
				u->Cell(it) = - u->Cell(it.Down());
				break;
			case boundaryRight:
				u->Cell(it.Left()) = 0.0;
				u->Cell(it) = 0.0;
				break;
			case boundaryLeft:
				u->Cell(it) = 0.0;
				break;
			case cornerTopRight:
				u->Cell(it.Left()) = 0.0;
				u->Cell(it) = - u->Cell(it.Down());
				break;
			case cornerTopLeft:
				u->Cell(it) = 0.0;
				break;
			case cornerBottomLeft:
				u->Cell(it) = 0.0;
				break;
			case cornerBottomRight:
				u->Cell(it.Left()) = 0.0;
				u->Cell(it) = -u->Cell(it.Top());
				break;
			}
			break;
		case 'O': //OUTFLOW
			switch (_type[it.Value()]) {
			case boundaryBottom:
				u->Cell(it) = u->Cell(it.Top());
				break;
			case boundaryTop:
				u->Cell(it) = u->Cell(it.Down());
				break;
			case boundaryRight:
				u->Cell(it.Left()) = u->Cell(it.Left().Left());
				u->Cell(it) = u->Cell(it.Left());
				break;
			case boundaryLeft:
				u->Cell(it) = u->Cell(it.Right());
				break;
			case cornerTopRight:
				u->Cell(it) = u->Cell(it.Down()) + u->Cell(it.Left());
				break;
			case cornerTopLeft:
				u->Cell(it) = u->Cell(it.Down()) + u->Cell(it.Right());
				break;
			case cornerBottomLeft:
				u->Cell(it) = u->Cell(it.Top()) + u->Cell(it.Right());
				break;
			case cornerBottomRight:
				u->Cell(it) = u->Cell(it.Top()) + u->Cell(it.Left());
				break;
			case inner:
				break;
			default:
				break;
			}
		case 'V': // Vertical INFLOW
			switch (_type[it.Value()]) {
			case boundaryBottom:
				u->Cell(it) = 2.0*_velocity[0] - u->Cell(it.Top());
				break;
			case boundaryTop:
				u->Cell(it) = 2.0*_velocity[0] - u->Cell(it.Down());
				break;
			case boundaryRight:
				u->Cell(it) = _value[it.Value()];
				u->Cell(it.Left()) = _value[it.Value()];
				break;
			case boundaryLeft:
				u->Cell(it) = _value[it.Value()];
				break;
			case cornerTopRight:
				if (_flag[it.Left().Value()] == ' ') u->Cell(it.Left()) = _value[it.Value()];
				u->Cell(it) = 2.0 * _velocity[1] - u->Cell(it.Down());
				break;
			case cornerTopLeft:
				u->Cell(it) = _value[it.Value()];
				break;
			case cornerBottomLeft:
				u->Cell(it) = _value[it.Value()];
				break;
			case cornerBottomRight:
				if (_flag[it.Left().Value()] == ' ') u->Cell(it.Left()) = _value[it.Value()];
				u->Cell(it) = 2.0 * _velocity[1] - u->Cell(it.Top());
				break;
			case inner:
				u->Cell(it) = _value[it.Value()];
				break;
			}
			break;
		case 'H':
			switch (_type[it.Value()]) {
			case boundaryBottom:
				u->Cell(it) = - u->Cell(it.Top());
				break;
			case boundaryTop:
				u->Cell(it) = -u->Cell(it.Down());
				break;
			case boundaryRight:
				u->Cell(it) = 0.0;
				u->Cell(it.Left()) = 0.0;
				break;
			case boundaryLeft:
				u->Cell(it) = 0.0;
				break;
			case cornerTopRight:
				u->Cell(it) = -u->Cell(it.Down());
				break;
			case cornerTopLeft:
				u->Cell(it.Left()) = 0.0;
				u->Cell(it) = - u->Cell(it.Down());
				break;
			case cornerBottomLeft:
				u->Cell(it.Left()) = 0.0;
				u->Cell(it) = -u->Cell(it.Top());
				break;
			case cornerBottomRight:
				u->Cell(it) = -u->Cell(it.Down());
				break;
			case inner:
				u->Cell(it) = 0.0;
				break;
			default:
				break;
			}
			break;
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
				v->Cell(it) = - v->Cell(it.Left());
				break;
			case boundaryBottom:
				v->Cell(it) = 0.0;
				break;
			case boundaryLeft:
				v->Cell(it) = - v->Cell(it.Right());
				break;
			case cornerTopRight:
				v->Cell(it.Down()) = 0.0;
				v->Cell(it) = - v->Cell(it.Left());
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
				v->Cell(it.Down()) = v->Cell(it.Down().Down());
				v->Cell(it) = v->Cell(it.Left());
				break;
			case cornerBottomRight:
				v->Cell(it) = v->Cell(it.Top());
				break;
			case cornerBottomLeft:
				v->Cell(it) = v->Cell(it.Top());
				break;
			case cornerTopLeft:
				v->Cell(it) = v->Cell(it.Right());
				v->Cell(it.Down()) = v->Cell(it.Down().Down());
				break;
			case inner:
				v->Cell(it) = 0.0;
				break;
			default:
				break;
			}
			break;
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
				v->Cell(it) = v->Cell(it.Down()) + v->Cell(it.Left());
				break;
			case cornerBottomRight:
				v->Cell(it) = v->Cell(it.Top()) + v->Cell(it.Left());
				break;
			case cornerBottomLeft:
				v->Cell(it) = v->Cell(it.Top()) + v->Cell(it.Right());
				break;
			case cornerTopLeft:
				v->Cell(it) = v->Cell(it.Down()) + v->Cell(it.Right());
				break;
			case inner:
				break;
			default:
				break;
			}
			break;
		case 'V': // Vertical INFLOW
			switch (_type[it.Value()]) {
			case boundaryTop:
				v->Cell(it) = 0.0;
				v->Cell(it.Down()) = 0.0;
				break;
			case boundaryRight:
				v->Cell(it) = - v->Cell(it.Left());
				break;
			case boundaryBottom:
				v->Cell(it) = 0.0;
				break;
			case boundaryLeft:
				v->Cell(it) = - v->Cell(it.Right());
				break;
			case cornerTopRight:
				v->Cell(it) = - v->Cell(it.Left());
				v->Cell(it.Down()) = 0.0;
				break;
			case cornerBottomRight:
				v->Cell(it) = - v->Cell(it.Left());
				break;
			case cornerBottomLeft:
				v->Cell(it) = - v->Cell(it.Right());
				break;
			case cornerTopLeft:
				v->Cell(it) = - v->Cell(it.Right());
				v->Cell(it.Down()) = 0.0;
				break;
			case inner:
				break;
			default:
				break;
			}
			break;
		case 'H': // Horizontal INFLOW
			switch (_type[it.Value()]) {
			case boundaryTop:
				v->Cell(it) = _value[it.Value()];
				v->Cell(it.Down()) = _value[it.Value()];
				break;
			case boundaryRight:
				v->Cell(it) = 2.0 * _velocity[1] - v->Cell(it.Left());
				break;
			case boundaryBottom:
				v->Cell(it) = _value[it.Value()];
				break;
			case boundaryLeft:
				v->Cell(it) = 2.0 * _velocity[1] - v->Cell(it.Right());
				break;
			case cornerTopRight:
				if(_flag[it.Down().Value()] == ' ') v->Cell(it.Down()) = _value[it.Value()]; 
				v->Cell(it) = 2.0 * _velocity[1] - v->Cell(it.Left());
				break;
			case cornerBottomRight:
				v->Cell(it) = _value[it.Value()];
				break;
			case cornerBottomLeft:
				v->Cell(it) = _value[it.Value()];
				break;
			case cornerTopLeft:
				if (_flag[it.Down().Value()] == ' ') v->Cell(it.Down()) = _value[it.Value()];
				v->Cell(it) = 2.0 * _velocity[1] - v->Cell(it.Right());
				break;
			case inner:
				break;
			default:
				break;
			}
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
				p->Cell(it) = 0.5 * (p->Cell(it.Left()) + p->Cell(it.Down()));
				break;
			case cornerBottomRight:
				p->Cell(it) = 0.5 * (p->Cell(it.Left()) + p->Cell(it.Top()));
				break;
			case cornerBottomLeft:
				p->Cell(it) = 0.5 * (p->Cell(it.Right()) + p->Cell(it.Top()));
				break;
			case cornerTopLeft:
				p->Cell(it) = 0.5 * (p->Cell(it.Right()) + p->Cell(it.Down()));
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
				p->Cell(it) = 2.0 * _pressure - p->Cell(it.Left());
				break;
			case boundaryBottom:
				p->Cell(it) = p->Cell(it.Top());
				break;
			case boundaryLeft:
				p->Cell(it) = 2.0 * _pressure - p->Cell(it.Right());
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
				p->Cell(it) = 2.0 * _pressure - p->Cell(it.Down());
				break;
			case boundaryRight:
				p->Cell(it) = p->Cell(it.Left());
				break;
			case boundaryBottom:
				p->Cell(it) = 2.0 * _pressure - p->Cell(it.Top());
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
			switch (_type[it.Value()]) {
			case boundaryTop:
				p->Cell(it) = - p->Cell(it.Down());
				break;
			case boundaryRight:
				p->Cell(it) = - p->Cell(it.Left());
				break;
			case boundaryBottom:
				p->Cell(it) = - p->Cell(it.Top());
				break;
			case boundaryLeft:
				p->Cell(it) = p->Cell(it.Right());
				break;
			case cornerTopRight:
				if (_flag[it.Down().Value()] == ' ') p->Cell(it) = p->Cell(it.Down());
				else p->Cell(it) = p->Cell(it.Left());
				break;
			case cornerBottomRight:
				if (_flag[it.Top().Value()] == ' ') p->Cell(it) = p->Cell(it.Top());
				else p->Cell(it) = p->Cell(it.Left());
				break;
			case cornerBottomLeft:
				if (_flag[it.Top().Value()] == ' ') p->Cell(it) = p->Cell(it.Top());
				else p->Cell(it) = p->Cell(it.Right());
				break;
			case cornerTopLeft:
				if (_flag[it.Down().Value()] == ' ') p->Cell(it) = p->Cell(it.Down());
				else p->Cell(it) = p->Cell(it.Right());
				break;
			case inner:
				p->Cell(it) = 0.0;
				break;
			default:
			break;
		}
		case 'V': // Vertical INFLOW
			break;
		case 'H': // Horizontal INFLOW
			break;
		default:
			break;
		}
	}
}
