#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include <cstdio>
#include <cstring>
#include <cstdlib>

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
    _size = multi_index_t(128, 128);
    _length  = multi_real_t (1.0, 1.0);
    _h  = multi_real_t (_length[0]/double(_size[0]), _length[1]/double(_size[1]));
    _velocity = multi_real_t (0.0 , 0.0);
    _pressure = 0.0;
    
	/*_size[0] = 128;
    _size[1] = 128;
	_length[0] = 1.0;
    _length[1] = 1.0;
	_h[0] = _length[0]/real_t(_size[0]);
	_h[1] = _length[1]/real_t(_size[1]);

	_velocity[0] = 0.0;
    _velocity[1] = 0.0;
	_pressure = 0.0;
    */
    
    
    Load("default.geom");
    
    //Load("actual.geom");
    
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
                _size[0] = inval_index[0];
                _size[1] = inval_index[1];
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


/// Updates the velocity field u
void Geometry::Update_U(Grid * u) const {
	BoundaryIterator bit = BoundaryIterator(new Geometry());
	
    bit.SetBoundary(0);
    // ecke oben links
    bit.First();
    u->Cell(bit.Left()) = 2.0 ;//- u->Cell(((bit.First()).Left()).Down())
	for(bit.First(); bit.Valid(); bit.Next()) {
		u->Cell(bit) = 2.0 - u->Cell(bit.Down());
	}
	
    bit.SetBoundary(1);
	for(bit.First(); bit.Valid(); bit.Next()) {
		u->Cell(bit.Left()) = 0.0;
	}
	
	bit.SetBoundary(2);
	for(bit.First(); bit.Valid(); bit.Next()) {
		u->Cell(bit) =  - u->Cell(bit.Top());
	}
	
	bit.SetBoundary(3);
    //ecke unten links
    bit.First();
    u->Cell(bit.Down()) = 0.0;
	for(bit.First(); bit.Valid(); bit.Next()) {
		u->Cell(bit) = 0.0;
	}
	
}

/// Updates the velocity field v
void Geometry::Update_V(Grid * v) const {
	BoundaryIterator bit = BoundaryIterator(new Geometry());
	bit.SetBoundary(0);
    
	for(bit.First(); bit.Valid(); bit.Next()) {
		v->Cell(bit.Down()) = 0.0;
	}
	bit.SetBoundary(1);
	for(bit.First(); bit.Valid(); bit.Next()) {
		v->Cell(bit) = - v->Cell(bit.Left());
	}
	bit.SetBoundary(2);
    //untere rechte ecke
    bit.First();
    v->Cell(bit.Right()) = 0.0;
	for(bit.First(); bit.Valid(); bit.Next()) {
		v->Cell(bit) = 0.0;
	}
	bit.SetBoundary(3);
    //untere linke ecke
    bit.First();
    v->Cell(bit.Down()) = 0.0;
	for(bit.First(); bit.Valid(); bit.Next()) {
		v->Cell(bit) = -v->Cell(bit.Right());
	}
}

/// Updates the pressure field p
void Geometry::Update_P(Grid * p) const {
    BoundaryIterator bit = BoundaryIterator(new Geometry());
	bit.SetBoundary(0);
    //hier werden die boundary conditions gesetzt. grad(p) soll null sein an allen rändern -> p_{0,i} = p{1,i} ...
	for(bit.First(); bit.Valid(); bit.Next()) {
		p->Cell(bit) = p->Cell(bit.Down());
	}
	bit.SetBoundary(1);
	for(bit.First(); bit.Valid(); bit.Next()) {
		p->Cell(bit) = p->Cell(bit.Left());
	}
	bit.SetBoundary(2);
	for(bit.First(); bit.Valid(); bit.Next()) {
		p->Cell(bit) = p->Cell(bit.Top());
	}
	bit.SetBoundary(3);
	for(bit.First(); bit.Valid(); bit.Next()) {
		p->Cell(bit) = p->Cell(bit.Right());
	}
}
