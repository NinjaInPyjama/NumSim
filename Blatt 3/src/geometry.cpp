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
        //funktioniert nur wenn das break drin ist sonst hängt er sich an der |#- geometry definition in default.geom auf.
        // |#- geometry definition in default.geom wird nicht eingelesen.
        if (strcmp(name,"geometry") == 0) {
            break;
        }
        
    }
	fclose(handle);
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


// /// Updates the velocity field u
// void Geometry::Update_U(Grid * u) const {
// 	// see script, p. 17
// 	BoundaryIterator bit = BoundaryIterator(this);
//     
//     // Iteration over right boundary
//     bit.SetBoundary(1);
// 	for(bit.First(); bit.Valid(); bit.Next()) {
// 		u->Cell(bit) = u->Cell(bit.Left());
// 	}
//     
//     // Iteration over left boundary
// 	bit.SetBoundary(3);
//     for(bit.First(); bit.Valid(); bit.Next()) {
// 		u->Cell(bit) = u->Cell(bit.Right());
// 	}
//     
//     // Iteration over upper boundary
//     bit.SetBoundary(0);
// 	for(bit.First(); bit.Valid(); bit.Next()) {
// 		u->Cell(bit) = - u->Cell(bit.Down()); //2.0
// 	}
// 	
//     // Iteration over lower boundary
// 	bit.SetBoundary(2);
// 	for(bit.First(); bit.Valid(); bit.Next()) {
// 		u->Cell(bit) = - u->Cell(bit.Top());
// 	}
// 	
// }
// 
// 
// 
// /// Updates the velocity field v
// void Geometry::Update_V(Grid * v) const {
// 	// see script, p. 17	
// 	BoundaryIterator bit = BoundaryIterator(this);
//     
//     // Iteration over right boundary
// 	bit.SetBoundary(1);
// 	for(bit.First(); bit.Valid(); bit.Next()) {
// 		v->Cell(bit) = v->Cell(bit.Left());
// 	}
// 	
// 	// Iteration over left boundary
// 	bit.SetBoundary(3);
// 	for(bit.First(); bit.Valid(); bit.Next()) {
// 		v->Cell(bit) = - v->Cell(bit.Right());
// 	}
//     
//     // Iteration over upper boundary
// 	bit.SetBoundary(0);
// 	for(bit.First(); bit.Valid(); bit.Next()) {
// 		
//         v->Cell(bit) = 0.0;
//         v->Cell(bit.Down()) = 0.0;
// 	}
//     
//     // Iteration over lower boundary
// 	bit.SetBoundary(2); // Lower right corner
// 	for(bit.First(); bit.Valid(); bit.Next()) {
// 		v->Cell(bit) = 0.0;
// 	}
//     
//     
//     
// 	
// }
// 
// /// Updates the pressure field p
// void Geometry::Update_P(Grid * p) const {
// 	// see script, p. 20 (p_{0,j} = p{1,j}), ...
// 	BoundaryIterator bit = BoundaryIterator(this);
//     
//     // Iteration over right boundary
// 	bit.SetBoundary(1);
// 	for(bit.First(); bit.Valid(); bit.Next()) {
// 		p->Cell(bit) = - p->Cell(bit.Left());
//     }
//     // Iteration over left boundary
// 	bit.SetBoundary(3);
// 	for(bit.First(); bit.Valid(); bit.Next()) {
// 		p->Cell(bit) = 2*_pressure - p->Cell(bit.Right());
// 	}
// 	
//     // Iteration over lower boundary
// 	bit.SetBoundary(2);
// 	for(bit.First(); bit.Valid(); bit.Next()) {
// 		p->Cell(bit) = p->Cell(bit.Top());
//     }
//     
// 	// Iteration over upper boundary
//     bit.SetBoundary(0);
//     for(bit.First(); bit.Valid(); bit.Next()) {
// 		p->Cell(bit) = p->Cell(bit.Down());
//     }
//     
//     
//     
// }



/// Updates the velocity field u
void Geometry::Update_U(Grid * u) const {
	// see script, p. 17
	BoundaryIterator bit = BoundaryIterator(this);
    
    // Iteration over right boundary
    bit.SetBoundary(1);
	for(bit.First(); bit.Valid(); bit.Next()) {
		u->Cell(bit.Left()) = 0.0;
	}
    
    // Iteration over left boundary
	bit.SetBoundary(3);
    bit.First();
    u->Cell(bit.Down()) = 0.0; // Lower left corner
    for(bit.First(); bit.Valid(); bit.Next()) {
		u->Cell(bit) = 0.0;
	}
    
    // Iteration over upper boundary
    bit.SetBoundary(0);
    bit.First();
	u->Cell(bit.Left()) = 2.0; // Upper left corner
	for(bit.First(); bit.Valid(); bit.Next()) {
		u->Cell(bit) = 2.0 - u->Cell(bit.Down()); //2.0
	}
	
    // Iteration over lower boundary
	bit.SetBoundary(2);
	for(bit.First(); bit.Valid(); bit.Next()) {
		u->Cell(bit) =  - u->Cell(bit.Top());
	}
	
}

/// Updates the velocity field v
void Geometry::Update_V(Grid * v) const {
	// see script, p. 17	
	BoundaryIterator bit = BoundaryIterator(this);
    
    // Iteration over upper boundary
	bit.SetBoundary(0);
	for(bit.First(); bit.Valid(); bit.Next()) {
		v->Cell(bit.Down()) = 0.0;
	}
    
    // Iteration over lower boundary
	bit.SetBoundary(2);
    bit.First();
    v->Cell(bit.Right()) = 0.0; // Lower right corner
	for(bit.First(); bit.Valid(); bit.Next()) {
		v->Cell(bit) = 0.0;
	}
    
    // Iteration over left boundary
	bit.SetBoundary(3);
    bit.First();
    v->Cell(bit.Down()) = 0.0; // Lower left corner
	for(bit.First(); bit.Valid(); bit.Next()) {
		v->Cell(bit) = - v->Cell(bit.Right());
	}
    
    // Iteration over right boundary
	bit.SetBoundary(1);
	for(bit.First(); bit.Valid(); bit.Next()) {
		v->Cell(bit) = - v->Cell(bit.Left());
	}
	
}

/// Updates the pressure field p
void Geometry::Update_P(Grid * p) const {
	// see script, p. 20 (p_{0,j} = p{1,j}), ...
	BoundaryIterator bit = BoundaryIterator(this);
    
	// Iteration over upper boundary
    bit.SetBoundary(0);
    for(bit.First(); bit.Valid(); bit.Next()) {
		p->Cell(bit) = p->Cell(bit.Down());
    }
    // Iteration over right boundary
	bit.SetBoundary(1);
	for(bit.First(); bit.Valid(); bit.Next()) {
		p->Cell(bit) = p->Cell(bit.Left());
    }
    // Iteration over lower boundary
	bit.SetBoundary(2);
	for(bit.First(); bit.Valid(); bit.Next()) {
		p->Cell(bit) = p->Cell(bit.Top());
    }
    
    // Iteration over left boundary
	bit.SetBoundary(3);
	for(bit.First(); bit.Valid(); bit.Next()) {
		p->Cell(bit) = p->Cell(bit.Right());
	}
}
