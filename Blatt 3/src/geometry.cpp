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
	_value = new real_t[_size[1] * _size[0]];

	_h = multi_real_t(_length[0] / (_size[0] - 2), _length[1] / (_size[1] - 2));
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
                _size[0] = inval_index[0]+2;
                _size[1] = inval_index[1]+2;
            }
            continue;
        }
        
        if (strcmp(name,"length") == 0) {
            if (fscanf(handle," %lf %lf\n",&inval_real[0],&inval_real[1])) {
                _length[0] = inval_real[0];
                _length[1] = inval_real[1];
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
        if (strcmp(name, "geometry") == 0) {
			if (fscanf(handle, " %s\n", name)) {
				if (strcmp(name, "free") == 0) {

					//char* _allflag = new char[_bsize[1]*_bsize[0]];
					_flag = new char[_size[1] * _size[0]];

					char* line = new char[_size[0]];
					for (int i = _size[1] - 1; i >= 0; i--) {

						fscanf(handle, "%[-HV| #IO]\n", line);
						//std::cout << line << std::endl;
						for (int j = 0; j < _size[0]; j++) {
							_flag[i*_size[0] + j] = line[j];
						}
					}
				}
			}
			continue;
		}
    }
	fclose(handle);
}

void Geometry::InitializeValues() {
	for (int i = 0; i < 4; i++) {
		index_t start_idx = -1;
		index_t firstID = _size[0] * (_size[1] - 1);
		index_t lastID = _size[0] * _size[1] - 1;
		index_t stepID = 1;
		switch (i) {
		case 1: 
			firstID = _size[0] * (_size[1] - 1);
			lastID = _size[0] * _size[1] - 1;
			stepID = 1;
			break;
		case 2:
			firstID = _size[0] * (_size[1] - 1);
			lastID = _size[0] * _size[1] - 1;
			stepID = 1;
			break;
		case 3:
			firstID = _size[0] * (_size[1] - 1);
			lastID = _size[0] * _size[1] - 1;
			stepID = 1;
			break;
		default:
			break;
		}
		for (int j = firstID; j <= lastID; j += stepID) {
			if (start_idx == -1 && (_flag[j] == 'H' || _flag[j] == 'V')) {
				start_idx = j;
			}
			else if (start_idx != -1 && (_flag[j] != 'H' && _flag[j] != 'V')) {
				index_t end_idx = j - stepID;

				//-4.0/(end_idx-start_idx+1.0)/(end_idx-start_idx+1.0)*velocity[1]*(y-start_idx+1.0/2.0)*(y-end_idx-1.0/2.0)
				
				for (int k = firstID; k <= lastID; k += stepID) {
					if (k <= end_idx && k >= start_idx) {
						_value[k] = -4.0 / ((end_idx - start_idx + 1.0)*(end_idx - start_idx + 1.0))*_velocity[(i + 1) % 2] * (real_t(k) - real_t(start_idx) + 1.0 / 2.0)*(real_t(k) - real_t(end_idx) - 1.0 / 2.0);
					}
				}
				start_idx = -1;
				end_idx = -1;
			}
		}
	}
}

/// Returns the number of cells in each dimension
const multi_index_t & Geometry::Size() const {
	return _size;
}

/// Returns the length of the domain
const multi_real_t & Geometry::Length() const {
	return _length;
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

/// Returns the flag array
const char * Geometry::Flag() const {
	return _flag;
}

/// Updates the velocity field u
void Geometry::Update_U(Grid * u) const {
	// see script, p. 17
	BoundaryIterator bit = BoundaryIterator(this);

	for (bit.First(); bit.Valid(); bit.Next()) {
		u->Cell(bit) = 0.0;
		switch (_flag[bit.Value()]) {
		case '#': // NOSLIP
			if (_flag[bit.Left().Value()] == ' ') { u->Cell(bit.Left()) = 0.0; u->Cell(bit) = 0.0; }
			if (_flag[bit.Top().Value()] == ' ') u->Cell(bit) = -u->Cell(bit.Top());
			if (_flag[bit.Down().Value()] == ' ') u->Cell(bit) = -u->Cell(bit.Down());
			if (_flag[bit.Right().Value()] == ' ') u->Cell(bit) = 0.0;
			break;
		case '-': // Horizontal SLIP
			if (_flag[bit.Top().Value()] == ' ') u->Cell(bit) = u->Cell(bit.Top());
			if (_flag[bit.Down().Value()] == ' ') u->Cell(bit) = u->Cell(bit.Down());
			if (_flag[bit.Left().Value()] == ' ') u->Cell(bit.Left()) = 0.0;
			if (_flag[bit.Right().Value()] == ' ') u->Cell(bit) = 0.0;
			break;
		case '|': // Vertical SLIP
			if (_flag[bit.Top().Value()] == ' ') u->Cell(bit) = -u->Cell(bit.Top());
			if (_flag[bit.Down().Value()] == ' ') u->Cell(bit) = -u->Cell(bit.Down());
			if (_flag[bit.Left().Value()] == ' ') u->Cell(bit.Left()) = 0.0;
			if (_flag[bit.Right().Value()] == ' ') u->Cell(bit) = 0.0;
			break;
		case 'O': // OUTFLOW
			if (_flag[bit.Top().Value()] == ' ') u->Cell(bit) = u->Cell(bit.Top());
			else if (_flag[bit.Down().Value()] == ' ') u->Cell(bit) = u->Cell(bit.Down());
			else if (_flag[bit.Left().Value()] == ' ') u->Cell(bit) = u->Cell(bit.Left());
			else if (_flag[bit.Right().Value()] == ' ') u->Cell(bit) = u->Cell(bit.Right());
			break;
		case 'V': // Vertical INFLOW
			if (_flag[bit.Top().Value()] == ' ') u->Cell(bit) = 2.0 * _velocity[0] - u->Cell(bit.Top());
			else if (_flag[bit.Down().Value()] == ' ') u->Cell(bit) = 2.0 * _velocity[0] - u->Cell(bit.Down());
			else if (_flag[bit.Left().Value()] == ' ') u->Cell(bit) = _value[bit.Value()];
			else if (_flag[bit.Right().Value()] == ' ') u->Cell(bit) = _value[bit.Value()];
			break;
		case 'H': // Horizontal INFLOW
			if (_flag[bit.Top().Value()] == ' ') u->Cell(bit) = -u->Cell(bit.Top());
			else if (_flag[bit.Down().Value()] == ' ') u->Cell(bit) = -u->Cell(bit.Down());
			else if (_flag[bit.Left().Value()] == ' ') u->Cell(bit.Left()) = 0.0;
			else if (_flag[bit.Right().Value()] == ' ') u->Cell(bit) = 0.0;
			break;
		default:
			break;
		}
	}
}



/// Updates the velocity field v
void Geometry::Update_V(Grid * v) const {
	// see script, p. 17	
	BoundaryIterator bit = BoundaryIterator(this);

	for (bit.First(); bit.Valid(); bit.Next()) {
		v->Cell(bit) = 0.0;
		switch (_flag[bit.Value()]) {
		case '#': // NOSLIP
			if (_flag[bit.Down().Value()] == ' ') { v->Cell(bit.Down()) = 0.0; v->Cell(bit) = 0.0; }
			if (_flag[bit.Left().Value()] == ' ') v->Cell(bit) = -v->Cell(bit.Left());
			if (_flag[bit.Right().Value()] == ' ') v->Cell(bit) = -v->Cell(bit.Right());
			if (_flag[bit.Top().Value()] == ' ') v->Cell(bit) = 0.0;
			break;
		case '-': // Horizontal SLIP
			if (_flag[bit.Left().Value()] == ' ') v->Cell(bit) = -v->Cell(bit.Left());
			if (_flag[bit.Right().Value()] == ' ') v->Cell(bit) = -v->Cell(bit.Right());
			if (_flag[bit.Down().Value()] == ' ') v->Cell(bit.Down()) = 0.0;
			if (_flag[bit.Top().Value()] == ' ') v->Cell(bit) = 0.0;
			break;
		case '|': // Vertical SLIP
			if (_flag[bit.Left().Value()] == ' ') v->Cell(bit) = v->Cell(bit.Left());
			if (_flag[bit.Right().Value()] == ' ') v->Cell(bit) = v->Cell(bit.Right());
			if (_flag[bit.Down().Value()] == ' ') v->Cell(bit.Down()) = 0.0;
			if (_flag[bit.Top().Value()] == ' ') v->Cell(bit) = 0.0;
			break;
		case 'O': // OUTFLOW
			if (_flag[bit.Left().Value()] == ' ') v->Cell(bit) = v->Cell(bit.Left());
			else if (_flag[bit.Right().Value()] == ' ') v->Cell(bit) = v->Cell(bit.Right());
			else if (_flag[bit.Down().Value()] == ' ') v->Cell(bit) = v->Cell(bit.Down());
			else if (_flag[bit.Top().Value()] == ' ') v->Cell(bit) = v->Cell(bit.Top());
			break;
		case 'V': // Vertical INFLOW
			if (_flag[bit.Left().Value()] == ' ') v->Cell(bit) = -v->Cell(bit.Left());
			else if (_flag[bit.Right().Value()] == ' ') v->Cell(bit) = -v->Cell(bit.Right());
			else if (_flag[bit.Down().Value()] == ' ') v->Cell(bit.Down()) = 0.0;
			else if (_flag[bit.Top().Value()] == ' ') v->Cell(bit) = 0.0;
			break;
		case 'H': // Horizontal INFLOW
			if (_flag[bit.Left().Value()] == ' ') v->Cell(bit) = 2.0 * _velocity[1] - v->Cell(bit.Left());
			else if (_flag[bit.Right().Value()] == ' ') v->Cell(bit) = 2.0 * _velocity[1] - v->Cell(bit.Right());
			else if (_flag[bit.Down().Value()] == ' ') v->Cell(bit) = _value[bit.Value()];
			else if (_flag[bit.Top().Value()] == ' ') v->Cell(bit) = _value[bit.Value()];
			break;
		default:
			break;
		}
	}
}

/// Updates the pressure field p
void Geometry::Update_P(Grid * p) const {
	// see script, p. 20 (p_{0,j} = p{1,j}), ...
	BoundaryIterator bit = BoundaryIterator(this);

	for (bit.First(); bit.Valid(); bit.Next()) {
		p->Cell(bit) = 0.0;
		switch (_flag[bit.Value()]) {
		case '#': // NOSLIP
			if (_flag[bit.Top().Value()] == ' ') p->Cell(bit) = p->Cell(bit.Top());
			if (_flag[bit.Down().Value()] == ' ') {
				if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Down());
				else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Down()));
			}
			if (_flag[bit.Left().Value()] == ' ') {
				if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Left());
				else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Left()));
			}
			if (_flag[bit.Right().Value()] == ' ') {
				if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Right());
				else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Right()));
			}
			break;
		case '-': // Horizontal SLIP
			if (_flag[bit.Top().Value()] == '-' || _flag[bit.Down().Value()] == '-') {
				if (_flag[bit.Right().Value()] == ' ') p->Cell(bit) = 2*_pressure - p->Cell(bit.Right());
				else p->Cell(bit) = 2.0 * _pressure - p->Cell(bit.Left());
			}
			else {
				if (_flag[bit.Top().Value()] == ' ') p->Cell(bit) = p->Cell(bit.Top());
				if (_flag[bit.Down().Value()] == ' ') {
					if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Down());
					else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Down()));
				}
				if (_flag[bit.Left().Value()] == ' ') {
					if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Left());
					else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Left()));
				}
				if (_flag[bit.Right().Value()] == ' ') {
					if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Right());
					else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Right()));
				}
			}
			break;
		case '|': // Vertical SLIP
			if (_flag[bit.Right().Value()] == '|' || _flag[bit.Left().Value()] == '|') {
				if (_flag[bit.Top().Value()] == ' ') p->Cell(bit) = 2 * _pressure - p->Cell(bit.Top());
				else p->Cell(bit) = 2.0 * _pressure - p->Cell(bit.Down());
			}
			else {
				if (_flag[bit.Top().Value()] == ' ') p->Cell(bit) = p->Cell(bit.Top());
				if (_flag[bit.Down().Value()] == ' ') {
					if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Down());
					else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Down()));
				}
				if (_flag[bit.Left().Value()] == ' ') {
					if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Left());
					else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Left()));
				}
				if (_flag[bit.Right().Value()] == ' ') {
					if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Right());
					else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Right()));
				}
			}
			break;
		case 'O': // OUTFLOW
			break;
		case 'V': // Vertical INFLOW
			if (_flag[bit.Top().Value()] == ' ') p->Cell(bit) = p->Cell(bit.Top());
			if (_flag[bit.Down().Value()] == ' ') {
				if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Down());
				else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Down()));
			}
			if (_flag[bit.Left().Value()] == ' ') {
				if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Left());
				else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Left()));
			}
			if (_flag[bit.Right().Value()] == ' ') {
				if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Right());
				else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Right()));
			}
			break;
		case 'H': // Horizontal INFLOW
			if (_flag[bit.Top().Value()] == ' ') p->Cell(bit) = p->Cell(bit.Top());
			if (_flag[bit.Down().Value()] == ' ') {
				if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Down());
				else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Down()));
			}
			if (_flag[bit.Left().Value()] == ' ') {
				if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Left());
				else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Left()));
			}
			if (_flag[bit.Right().Value()] == ' ') {
				if (p->Cell(bit) == 0.0) p->Cell(bit) = p->Cell(bit.Right());
				else p->Cell(bit) = 0.5*(p->Cell(bit) + p->Cell(bit.Right()));
			}
			break;
		default:
			break;
		}
	}
}

/*
/// Updates the velocity field u
void Geometry::Update_U(Grid * u) const {
// see script, p. 17
BoundaryIterator bit = BoundaryIterator(this);

// Iteration over right boundary
bit.SetBoundary(1);
for(bit.First(); bit.Valid(); bit.Next()) {
u->Cell(bit) = u->Cell(bit.Left());
}

// Iteration over left boundary
bit.SetBoundary(3);
for(bit.First(); bit.Valid(); bit.Next()) {
u->Cell(bit) = u->Cell(bit.Right());
}

// Iteration over upper boundary
bit.SetBoundary(0);
for(bit.First(); bit.Valid(); bit.Next()) {
u->Cell(bit) = - u->Cell(bit.Down()); //2.0
}

// Iteration over lower boundary
bit.SetBoundary(2);
for(bit.First(); bit.Valid(); bit.Next()) {
u->Cell(bit) = - u->Cell(bit.Top());
}

}



/// Updates the velocity field v
void Geometry::Update_V(Grid * v) const {
// see script, p. 17
BoundaryIterator bit = BoundaryIterator(this);

// Iteration over right boundary
bit.SetBoundary(1);
for(bit.First(); bit.Valid(); bit.Next()) {
v->Cell(bit) = v->Cell(bit.Left());
}

// Iteration over left boundary
bit.SetBoundary(3);
for(bit.First(); bit.Valid(); bit.Next()) {
v->Cell(bit) = - v->Cell(bit.Right());
}

// Iteration over upper boundary
bit.SetBoundary(0);
for(bit.First(); bit.Valid(); bit.Next()) {

v->Cell(bit) = 0.0;
v->Cell(bit.Down()) = 0.0;
}

// Iteration over lower boundary
bit.SetBoundary(2); // Lower right corner
for(bit.First(); bit.Valid(); bit.Next()) {
v->Cell(bit) = 0.0;
}




}

/// Updates the pressure field p
void Geometry::Update_P(Grid * p) const {
// see script, p. 20 (p_{0,j} = p{1,j}), ...
BoundaryIterator bit = BoundaryIterator(this);

// Iteration over right boundary
bit.SetBoundary(1);
for(bit.First(); bit.Valid(); bit.Next()) {
p->Cell(bit) = - p->Cell(bit.Left());
}
// Iteration over left boundary
bit.SetBoundary(3);
for(bit.First(); bit.Valid(); bit.Next()) {
p->Cell(bit) = 2*_pressure - p->Cell(bit.Right());
}

// Iteration over lower boundary
bit.SetBoundary(2);
for(bit.First(); bit.Valid(); bit.Next()) {
p->Cell(bit) = p->Cell(bit.Top());
}

// Iteration over upper boundary
bit.SetBoundary(0);
for(bit.First(); bit.Valid(); bit.Next()) {
p->Cell(bit) = p->Cell(bit.Down());
}
}
*/