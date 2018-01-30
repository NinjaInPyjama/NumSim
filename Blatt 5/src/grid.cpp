#include "grid.hpp"
#include "iterator.hpp"

/// Constructs a grid based on a geometry
Grid::Grid(const Geometry * geom) {
	_data = new real_t[geom->Size()[0] * geom->Size()[1]];
	_offset = multi_real_t(0.0, 0.0);
	_geom = geom;
}

/// Constructs a grid based on a geometry with an offset
// @param geom   Geometry information
// @param offset distance of staggered grid point to cell's anchor point;
//               (anchor point = lower Left corner)
Grid::Grid(const Geometry * geom, const multi_real_t & offset) {
	_data = new real_t[geom->Size()[0] * geom->Size()[1]];
	_offset = offset;
	_geom = geom;
}

Grid::Grid() {}



/// Deletes the grid
Grid::~Grid() {}

/// Prints the values of grid
void Grid::print() const {
	std::cout.precision(2);
	for (int i = _geom->Size()[1] - 1; i >= 0; i--) {
		for (int j = 0; j < _geom->Size()[0]; j++) {
			std::cout << " " << _data[i*_geom->Size()[0] + j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

///  Initializes the grid with a value
void Grid::Initialize(const real_t & value) {
    index_t num_cells = _geom->Size()[0] * _geom->Size()[1];
    for (index_t i = 0; i < num_cells; i++) {
		_data[i] = value;
	}
}


/// Write access to the grid cell at position [it]
real_t & Grid::Cell(const Iterator & it) {
	return _data[it.Value()];
}

/// Read access to the grid cell at position [it]
const real_t & Grid::Cell(const Iterator & it) const {
	return _data[it.Value()];
}


/// Interpolate the value at an arbitrary position by bilinear interpolation
real_t Grid::Interpolate(const multi_real_t & pos) const {
//     multi_real_t newpos = multi_real_t(pos[0], pos[1]);
//     if(pos[0] > _geom->Length()[0]) newpos[0] = _geom->Length()[0];
//     if(pos[0] < 0.0) newpos[0] = 0.0;
//     if(pos[1] > _geom->Length()[1]) newpos[1] = _geom->Length()[1];
//     if(pos[1] < 0.0) newpos[1] = 0.0;
    
    real_t pos_x = pos[0] * (_geom->Size()[0] - 2) / _geom->Length()[0] + 1 - _offset[0];
    real_t pos_y = pos[1] * (_geom->Size()[1] - 2) / _geom->Length()[1] + 1 - _offset[1];
    
    if(pos_x >= _geom->Size()[0]-1) pos_x = _geom->Size()[0] - 1 - 0.00000001;
    if(pos_x < 0.0) pos_x = 0;
    if(pos_y >= _geom->Size()[1]-1) pos_y = _geom->Size()[1] - 1 - 0.00000001;
    if(pos_y < 0.0) pos_y = 0;
    
	index_t index_x = (index_t)pos_x;
	index_t index_y = (index_t)pos_y;

	// Lower left data point
    real_t val_ll = _data[index_x + index_y*_geom->Size()[0]];
	// Lower Right data point
	real_t val_lr = _data[index_x + 1 + index_y*_geom->Size()[0]];
	// Upper left data point
    real_t val_ul = _data[index_x + (index_y + 1)*_geom->Size()[0]];
	// Upper Right data point
    real_t val_ur = _data[index_x + 1 + (index_y + 1)*_geom->Size()[0]];
    
	// Proportion in x-dim
    real_t prop_x = pos_x - (real_t)index_x; 
	// Proportion in y-dim
    real_t prop_y = pos_y - (real_t)index_y;

	//std::cout << pos_x << " " << pos_y << " " << index_x << " " << index_y << " " << prop_x << " " << prop_y << std::endl;

	return (val_ll*(1.0 - prop_x) + prop_x*val_lr)*(1.0-prop_y) + prop_y*( val_ul*(1.0 - prop_x) + prop_x*val_ur );
}


/// Computes the left-sided difference quatient in x-dim at [it]
real_t Grid::dx_l(const Iterator & it) const {
	return (Cell(it) - Cell(it.Left()))/_geom->Mesh()[0];
}

/// Computes the Right-sided difference quatient in x-dim at [it]
real_t Grid::dx_r(const Iterator & it) const {
	return (Cell(it.Right()) - Cell(it)) / _geom->Mesh()[0];
}

/// Computes the left-sided difference quatient in y-dim at [it]
real_t Grid::dy_l(const Iterator & it) const {
	return (Cell(it) - Cell(it.Down())) / _geom->Mesh()[1];
}

/// Computes the Right-sided difference quatient in x-dim at [it]
real_t Grid::dy_r(const Iterator & it) const {
	return (Cell(it.Top()) - Cell(it)) / _geom->Mesh()[1];
}

/// Computes the central difference quatient of 1st order in x-dim at [it]
real_t Grid::dx_central(const Iterator & it) const {
	return (Cell(it.Right()) - Cell(it.Left())) / (2.0 * _geom->Mesh()[0]);
}

/// Computes the central difference quatient of 1st order in y-dim at [it]
real_t Grid::dy_central(const Iterator & it) const {
	return (Cell(it.Top()) - Cell(it.Down())) / (2.0 * _geom->Mesh()[1]);
}

/// Computes the central difference quatient of 2nd order in x-dim at [it]
real_t Grid::dxx(const Iterator & it) const {
	return (Cell(it.Right()) - 2.0*Cell(it) + Cell(it.Left())) / (_geom->Mesh()[0] * _geom->Mesh()[0]) ;
}

/// Computes the central difference quatient of 2nd order in y-dim at [it]
real_t Grid::dyy(const Iterator & it) const {
	return (Cell(it.Top()) - 2.0*Cell(it) + Cell(it.Down())) / (_geom->Mesh()[1] * _geom->Mesh()[1]);
}


/// Computes u*du/dx with the donor cell method
real_t Grid::DC_udu_x(const Iterator & it, const real_t & alpha) const {
	return 0;
}

/// Computes v*du/dy with the donor cell method
real_t Grid::DC_vdu_y(const Iterator & it, const real_t & alpha, const Grid * v) const {
	return 0;
}

/// Computes u*dv/dx with the donor cell method
real_t Grid::DC_udv_x(const Iterator & it, const real_t & alpha, const Grid * u) const {
	return 0;
}

/// Computes v*dv/dy with the donor cell method
real_t Grid::DC_vdv_y(const Iterator & it, const real_t & alpha) const {
	return 0;
}



/// Computes du^2/dx with the donor cell method
real_t Grid::DC_du2_x(const Iterator & it, const real_t & alpha) const {
	// see script, p.22

	const real_t dx = _geom->Mesh()[0];

	// Value of u at iterator cell (u_{i,j})
	const real_t val_u = Cell(it);
	// Value of u at the Right neighbor of the iterator cell (u_{i+1,j})
	const real_t val_u_r = Cell(it.Right());
	// Value of u at the left neighbor of the iterator cell (u_{i-1,j})
	const real_t val_u_l = Cell(it.Left());
	// Interpolated value of u between this and its right neighbor cell (u_{i+1/2,j})
	const real_t val_u_cr = (val_u_r + val_u) / 2.0;
	// Interpolated value of u between this and its left neighbor cell (u_{i-1/2,j})
	const real_t val_u_cl = (val_u + val_u_l) / 2.0;
	return (val_u_cr * val_u_cr - val_u_cl * val_u_cl) / dx
			+ alpha * (abs(val_u_cr) * (val_u - val_u_r) / 2.0 - abs(val_u_cl) * (val_u_l - val_u) / 2.0) / dx;
}

/// Computes dv^2/dy with the donor cell method
real_t Grid::DC_dv2_y(const Iterator & it, const real_t & alpha) const {
	// see script, p.22

	const real_t dy = _geom->Mesh()[1];

	// Value of v at iterator cell (v_{i,j})
	const real_t val_v = Cell(it);
	// Value of v at the upper neighbor (Top) of the iterator cell (v_{i,j+1})
	const real_t val_v_t = Cell(it.Top());
	// Value of v at the lower neighbor (Down) of the iterator cell (v_{i,j-1})
	const real_t val_v_d = Cell(it.Down());
	// Interpolated value of v between this and its upper neighbor cell (v_{i,j+1/2})
	const real_t val_v_ct = (val_v_t + val_v) / 2.0;
	// Interpolated value of v between this and its lower neighbor cell (v_{i,j-1/2})
	const real_t val_v_cd = (val_v + val_v_d) / 2.0;
	return (val_v_ct * val_v_ct - val_v_cd * val_v_cd) / dy
			+ alpha * (abs(val_v_ct) * (val_v - val_v_t) / 2.0 - abs(val_v_cd) * (val_v_d - val_v) / 2.0) / dy;
}

/// Computes d(uv)/dx with the donor cell method
real_t Grid::DC_duv_x(const Iterator & it, const real_t & alpha, const Grid * u) const {
	// see script, p.22

	const real_t dx = _geom->Mesh()[0];

	// Value of v at iterator cell (v_{i,j})
	const real_t val_v = Cell(it);
	// Value of v at the right neighbor of the iterator cell (v_{i+1,j})
	const real_t val_v_r = Cell(it.Right());
	// Value of u at the left neighbor of the iterator cell (v_{i-1,j})
	const real_t val_v_l = Cell(it.Left());

	// Interpolated value of u between this and its upper neighbor cell (u_{i,j+1/2})
	const real_t val_u_ct = (u->Cell(it.Top()) + u->Cell(it)) / 2.0;
	// Interpolated value of u between the iterator cell's left and its upper left neighbor cell (u_{i-1,j+1/2})
	const real_t val_u_ctl = (u->Cell(it.Left().Top()) + u->Cell(it.Left())) / 2.0;

	return (val_u_ct * (val_v + val_v_r) / 2.0 - val_u_ctl * (val_v_l + val_v) / 2.0) / dx +
			+ alpha * (abs(val_u_ct) * (val_v - val_v_r) / 2.0 - abs(val_u_ctl) * (val_v_l - val_v) / 2.0) / dx;
}

/// Computes d(uv)/dy with the donor cell method
real_t Grid::DC_duv_y(const Iterator & it, const real_t & alpha, const Grid * v) const {
	// see script, p.22

	const real_t dy = _geom->Mesh()[1];

	// Value of u at iterator cell (u_{i,j})
	const real_t val_u = Cell(it);
	// Value of u at the upper neighbor of the iterator cell (u_{i,j+1})
	const real_t val_u_t = Cell(it.Top());
	// Value of u at the lower neighbor of the iterator cell (u_{i,j-1})
	const real_t val_u_d = Cell(it.Down());

	// Interpolated value of v between this and its right neighbor cell (v_{i+1/2,j})
	const real_t val_v_cr = (v->Cell(it.Right()) + v->Cell(it)) / 2.0;
	// Interpolated value of v between the iterator cell's lower and its lower right neighbor cell (v_{i+1/2,j-1})
	const real_t val_v_cdr = (v->Cell(it.Right().Down()) + v->Cell(it.Down())) / 2.0;

	return (val_v_cr * (val_u + val_u_t) / 2.0 - val_v_cdr * (val_u_d + val_u) / 2.0) / dy 
			+ alpha * (abs(val_v_cr) * (val_u - val_u_t) / 2.0 - abs(val_v_cdr) * (val_u_d - val_u) / 2.0) / dy;
}


/// Returns the maximal value of the grid
real_t Grid::Max() const {
	real_t max = _data[0];
	const index_t num_cells = _geom->Size()[0] * _geom->Size()[1];
	for (index_t i = 1; i < num_cells; i++) {
		if(max < _data[i]) max = _data[i];
	}
	return max;
}

/// Returns the minimal value of the grid
real_t Grid::Min() const {
	real_t min = _data[0];
	const index_t num_cells = _geom->Size()[0] * _geom->Size()[1];
	for (index_t i = 1; i < num_cells; i++) {
		if (min > _data[i]) min = _data[i];
	}
	return min;
}

/// Returns the absolute maximal value
real_t Grid::AbsMax() const {
	real_t max = std::abs(_data[0]);
	const index_t num_cells = _geom->Size()[0] * _geom->Size()[1];
	for (index_t i = 1; i < num_cells; i++) {
		if (max < std::abs(_data[i])) max = std::abs(_data[i]);
	}
	return max;
}

/// Returns a pointer to the raw data
real_t * Grid::Data() {
	return _data;
}





MultiGrid::MultiGrid(const Geometry *geom) {
    _geom = geom;
    _data = new real_t[geom->Size()[0] * geom->Size()[1]];
    _cellSize = new index_t[geom->Size()[0] * geom->Size()[1]];
    Initialize(0.0);
    for(int i=0; i< (geom->Size()[0] * geom->Size()[1]); i++) _cellSize[i] = 1;
}

  
MultiGrid::MultiGrid(const Geometry *geom, index_t* cellSize, real_t* data ){
    _geom = geom;
    _data = data;
    _cellSize = cellSize;
}
  
MultiGrid::MultiGrid(Grid *grid, const Geometry *geom) {
    _geom = geom;
    _data = grid->Data();
    _data = new real_t[geom->Size()[0] * geom->Size()[1]];
    _cellSize = new index_t[geom->Size()[0] * geom->Size()[1]];
    
    for(int i=0; i<geom->Size()[0] * geom->Size()[1]; i++) {
        _cellSize[i] = 1;
    }
}


MultiGrid::MultiGrid() {}

  /// Deletes the multigrid
MultiGrid::~MultiGrid() {}

real_t &MultiGrid::Cell(const index_t index) {
    return _data[index];
}

void MultiGrid::printCellSize() const{
    for (int i = _geom->Size()[1] - 1; i >= 0; i--) {
        for (int j = 0; j < _geom->Size()[0]; j++) {
            std::cout << " " << _cellSize[i*_geom->Size()[0] + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

const index_t &MultiGrid::CellSize(const index_t index) const{
    return _cellSize[index];
}
  
  /// Write access to the grid cell at position [it]
//real_t &MultiGrid::Cell(const Iterator &it) {
    
//}
  /// Read access to the grid cell at position [it]
//const real_t &MultiGrid::Cell(const Iterator &it) const{
    
//}

real_t MultiGrid::dxx(const MGIterator &it){
    if(_cellSize[it.Value()] != 0) {
      real_t h = _cellSize[it.Value()]*_geom->Mesh()[0];
      if(_cellSize[it.MGLeft().Value()] == 0.5*h) {
        if(_cellSize[it.MGRight().Value()] == 0.5*h) {
          // Left and right finer
          return 16/9 * (0.5 * (Cell(it.MGLeft()) + Cell(it.MGLeft().MGTop())) - 2.0 * Cell(it) + 0.5 * (Cell(it.MGRight()) + Cell(it.MGRight().MGTop()))) / (h*h);
        }
        else {
          // Left finer, right same or coarser
          // from matlab dxx = 2 * (4*(4*vl - 7*vc + 3*vr))/(21*h^2)
          return 8/21 * (4.0 * 0.5 * (Cell(it.MGLeft()) + Cell(it.MGLeft().MGTop())) - 7.0 * Cell(it) + 3.0 * Cell(it.MGRight())) / (h*h);
        }
      }
      else {
        if(_cellSize[it.MGRight().Value()] == 0.5*h) {
          // Right same or coarser, left finer
          return 8/21 * (3.0 * Cell(it.MGLeft()) - 7.0 * Cell(it) + 4.0 * 0.5 * (Cell(it.MGRight()) + Cell(it.MGRight().MGTop()))) / (h*h);
        }
        else {
          // Left and right same or coarser
          return (Cell(it.MGRight()) - 2.0 * Cell(it) + Cell(it.MGLeft())) / (h*h);
        }
      }
    }
    else {
        
        std::cout << "access to zero field" << std::endl;
        return 0.0;
    }
}
  
real_t MultiGrid::dyy(const MGIterator &it){
    if(_cellSize[it.Value()] != 0) {
      real_t h = _cellSize[it.Value()]*_geom->Mesh()[1];
      if(_cellSize[it.MGTop().Value()] == 0.5*h) {
        if(_cellSize[it.MGDown().Value()] == 0.5*h) {
          // Top and bottom finer
          return 16/9 * (0.5 * (Cell(it.MGTop()) + Cell(it.MGTop().MGRight())) - 2.0 * Cell(it) + 0.5 * (Cell(it.MGDown()) + Cell(it.MGDown().MGRight()))) / (h*h);
        }
        else {
          // Top finer, bottom same or coarser
          // from matlab dxx = 2 * (4*(4*vt - 7*vc + 3*vb))/(21*h^2)
          return 8/21 * (4.0 * 0.5 * (Cell(it.MGTop()) + Cell(it.MGTop().MGRight())) - 7.0 * Cell(it) + 3.0 * Cell(it.MGDown())) / (h*h);
        }
      }
      else {
        if(_cellSize[it.MGDown().Value()] == 0.5*h) {
          // Top same or coarser, bottom finer
          return 8/21 * (3.0 * Cell(it.MGTop()) - 7.0 * Cell(it) + 4.0 * 0.5 * (Cell(it.MGDown()) + Cell(it.MGDown().MGRight()))) / (h*h);
        }
        else {
          // Top and bottom same or coarser
          return (Cell(it.MGTop()) - 2.0 * Cell(it) + Cell(it.MGDown())) / (h*h);
        }
      }
    }
    else {
        
        std::cout << "access to zero field" << std::endl;
        return 0.0;
    }
}
  
MultiGrid * MultiGrid::restrict(const index_t resSize) {
    index_t* newCellSize = new index_t[_geom->Size()[0]*_geom->Size()[1]];
    real_t* newData = new real_t[_geom->Size()[0]*_geom->Size()[1]];
    for(int i=0; i<_geom->Size()[0] * _geom->Size()[1]; i++) {
        newCellSize[i] = _cellSize[i];
        newData[i] = _data[i];
    }
    
    
    MGInteriorIterator it = MGInteriorIterator(_geom, new MultiGrid(_geom, newCellSize, newData), resSize);
    
    for(it.First(); it.Valid(); it.Next()) {
      if(_cellSize[it.MGTop().Value()] == resSize && _geom->Flag()[it.MGTop().Value()] == ' ' 
         && _cellSize[it.MGRight().Value()] == resSize & _geom->Flag()[it.MGRight().Value()] == ' ' 
         && _cellSize[it.MGTop().MGRight().Value()] == resSize && _geom->Flag()[it.MGTop().MGRight().Value()] == ' '
         && _cellSize[it.MGLeft().Value()] == resSize
         && _cellSize[it.MGLeft().MGTop().Value()] == resSize
         && _cellSize[it.MGDown().Value()] == resSize
         && _cellSize[it.MGDown().MGRight().Value()] == resSize
         && _cellSize[it.MGTop().MGTop().Value()] == resSize
         && _cellSize[it.MGTop().MGTop().MGRight().Value()] == resSize
         && _cellSize[it.MGRight().MGRight().Value()] == resSize
         && _cellSize[it.MGRight().MGRight().MGTop().Value()] == resSize) {
        
         Cell(it) =  0.25*(Cell(it) + Cell(it.MGTop()) + Cell(it.MGRight()) + Cell(it.MGTop().MGRight()));
         Cell(it.Value()+2*resSize-1) =  0.25*(Cell(it) + Cell(it.MGTop()) + Cell(it.MGRight()) + Cell(it.MGTop().MGRight()));
         Cell(it.Value()+2*resSize*(_geom->Size()[0]-1)) =  0.25*(Cell(it) + Cell(it.MGTop()) + Cell(it.MGRight()) + Cell(it.MGTop().MGRight()));
         // add here the other 4 points
         
         int index = it.Value();
         int indexTop = it.MGTop().Value();
         int indexRight = it.MGRight().Value();
         int indexTopRight = it.MGTop().MGRight().Value();
        
         newCellSize[index] = 2*resSize;
         newCellSize[indexTop] = 0;
         newCellSize[indexRight] = 0;
         newCellSize[indexTopRight] = 0;
      }
    }
    
    MGBoundaryIterator bit = MGBoundaryIterator(_geom, new MultiGrid(_geom, newCellSize, newData), resSize);

    bit.SetBoundary(0);
    for(bit.First(); bit.Valid(); bit.Next()) {
      if(_cellSize[bit.MGTop().Value()] == resSize & _geom->Flag()[bit.MGTop().Value()] == _geom->Flag()[bit.Value()]
         && _cellSize[bit.MGLeft().Value()] == resSize
         && _cellSize[bit.MGLeft().MGTop().Value()] == resSize) {
         
         //write2Cell(bit, 0.5*(Cell(bit) + Cell(bit.Top())));
        
         int index = bit.Value();
         int indexRight = bit.MGTop().Value();
        
         newCellSize[index] = 2*resSize;
         newCellSize[indexRight] = 0;
       }
    }
    
    bit.SetBoundary(1);
    for(bit.First(); bit.Valid(); bit.Next()) {
      if(_cellSize[bit.MGRight().Value()] == resSize & _geom->Flag()[bit.MGRight().Value()] == _geom->Flag()[bit.Value()]
         && _cellSize[bit.MGTop().Value()] == resSize
         && _cellSize[bit.MGTop().MGRight().Value()] == resSize) {
         
         //write2Cell(bit, 0.5*(Cell(bit) + Cell(bit.Right())));
         
         int index = bit.Value();
         int indexRight = bit.MGRight().Value();
        
         newCellSize[index] = 2*resSize;
         newCellSize[indexRight] = 0;
       }
    }
    
    bit.SetBoundary(2);
    for(bit.First(); bit.Valid(); bit.Next()) {
      if(_cellSize[bit.MGTop().Value()] == resSize & _geom->Flag()[bit.MGTop().Value()] == _geom->Flag()[bit.Value()]
         && _cellSize[bit.MGRight().Value()] == resSize
         && _cellSize[bit.MGRight().MGTop().Value()] == resSize) {
         
         //write2Cell(bit, 0.5*(Cell(bit) + Cell(bit.Top())));
        
         int index = bit.Value();
         int indexRight = bit.MGTop().Value();
        
         newCellSize[index] = 2*resSize;
         newCellSize[indexRight] = 0;
       }
    }
    
    bit.SetBoundary(3);
    for(bit.First(); bit.Valid(); bit.Next()) {
      if(_cellSize[bit.MGRight().Value()] == resSize & _geom->Flag()[bit.MGRight().Value()] == _geom->Flag()[bit.Value()]
         && _cellSize[bit.MGDown().Value()] == resSize
         && _cellSize[bit.MGDown().MGRight().Value()] == resSize) {
         
         //write2Cell(bit, 0.5*(Cell(bit) + Cell(bit.Right())));
        
         int index = bit.Value();
         int indexRight = bit.MGRight().Value();
        
         newCellSize[index] = 2*resSize;
         newCellSize[indexRight] = 0;
       }
    }
    _geom->Update_P(new MultiGrid(_geom, newCellSize, newData));
    
    return new MultiGrid(_geom, newCellSize, newData);
}
  
MultiGrid * MultiGrid::prolong(const index_t intSize) {
    index_t* newCellSize = new index_t[_geom->Size()[0]*_geom->Size()[1]];
    real_t* newData = new real_t[_geom->Size()[0]*_geom->Size()[1]];
    for(int i=0; i<_geom->Size()[0] * _geom->Size()[1]; i++) {
        newCellSize[i] = _cellSize[i];
        newData[i] = _data[i];
    }
    
    
    MGInteriorIterator it = MGInteriorIterator(_geom, new MultiGrid(_geom, newCellSize, newData), intSize);
    for(it.First(); it.Valid(); it.Next()) {
       if(_cellSize[it.Value()] == intSize && _geom->Flag()[it.Value()] == ' ') {
         
         newCellSize[it.Value()] = intSize/2;
         newCellSize[it.Value()+intSize/2] = intSize/2;
         newCellSize[it.Value()+intSize/2*_geom->Size()[0]] = intSize/2;
         newCellSize[it.Value()+intSize/2*(1+_geom->Size()[0])] = intSize/2;
         
         Cell(it.MGTop()) =  Cell(it);
         Cell(it.MGRight()) = Cell(it);
         Cell(it.MGTop().MGRight()) = Cell(it);
         
         Cell(it.Value()+intSize/2-1) = Cell(it);
         Cell(it.MGTop().Value()+intSize/2-1) = Cell(it);
         Cell(it.MGRight().Value()+intSize/2-1) = Cell(it);
         Cell(it.MGTop().MGRight().Value()+intSize/2-1) = Cell(it);
         
         Cell(it.Value()+intSize/2*(_geom->Size()[0]-1)) = Cell(it);
         Cell(it.MGTop().Value()+intSize/2*(_geom->Size()[0]-1)) = Cell(it);
         Cell(it.MGRight().Value()+intSize/2*(_geom->Size()[0]-1)) = Cell(it);
         Cell(it.MGTop().MGRight().Value()+intSize/2*(_geom->Size()[0]-1)) = Cell(it);
       }
    }
    
    MGBoundaryIterator bit = MGBoundaryIterator(_geom, new MultiGrid(_geom, newCellSize, newData), intSize);

    bit.SetBoundary(0);
    for(bit.First(); bit.Valid(); bit.Next()) {
      if(_cellSize[bit.Value()] == intSize) {
        
        newCellSize[bit.Value()] = intSize/2;
        newCellSize[bit.Value() + intSize/2*_geom->Size()[0]] = intSize/2;
      }
    }
    
    bit.SetBoundary(1);
    for(bit.First(); bit.Valid(); bit.Next()) {
      if(_cellSize[bit.Value()] == intSize) {
        
        newCellSize[bit.Value()] = intSize/2;
        newCellSize[bit.Value() + intSize/2] = intSize/2;
      }
    }
    
    bit.SetBoundary(2);
    for(bit.First(); bit.Valid(); bit.Next()) {
      if(_cellSize[bit.Value()] == intSize) {
        
        newCellSize[bit.Value()] = intSize/2;
        newCellSize[bit.Value() + intSize/2*_geom->Size()[0]] = intSize/2;
      }
    }
    
    bit.SetBoundary(3);
    for(bit.First(); bit.Valid(); bit.Next()) {
      if(_cellSize[bit.Value()] == intSize) {
        
        newCellSize[bit.Value()] = intSize/2;
        newCellSize[bit.Value() + intSize/2] = intSize/2;
      }
    }
    
    return new MultiGrid(_geom, newCellSize, newData);
  }

