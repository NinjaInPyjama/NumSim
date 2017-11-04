#include "grid.hpp"

/// Constructs a grid based on a geometry
Grid::Grid(const Geometry * geom) {
	_data = new real_t[geom->Size()[0] * geom->Size()[1]];
	_offset = multi_real_t(0.0, 0.0);
	_geom = geom;
}

/// Constructs a grid based on a geometry with an offset
// @param geom   Geometry information
// @param offset distance of staggered grid point to cell's anchor point;
//               (anchor point = lower left corner)
Grid::Grid(const Geometry * geom, const multi_real_t & offset) {
	_data = new real_t[geom->Size()[0] * geom->Size()[1]];
	_offset = offset;
	_geom = geom;
}


/// Deletes the grid
Grid::~Grid() {}


///     Initializes the grid with a value
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


/// Interpolate the value at an arbitrary position
/// bilinear interpolation, 
real_t Grid::Interpolate(const multi_real_t & pos) const {
    index_t pos_x = (index_t)(pos[0] - _offset[0]); //integer value in x direction of lower left corner
    index_t pos_y = (index_t)(pos[1] - _offset[1]); //integer value in y direction of lower left corner
    
	// Lower left data point
    real_t val_ll = _data[pos_x + pos_y*_geom->Size()[0]];
	// Lower right data point
	real_t val_lr = _data[pos_x + 1 + pos_y*_geom->Size()[0]];
	// Upper left data point
    real_t val_ul = _data[pos_x + (pos_y + 1)*_geom->Size()[0]];
	// Upper right data point
    real_t val_ur = _data[pos_x + 1 + (pos_y + 1)*_geom->Size()[0]];
    
	// Proportion in x-dim
    real_t prop_x = pos[0] - _offset[0] - (real_t)pos_x; 
	// Proportion in y-dim
    real_t prop_y = pos[1] - _offset[1] - (real_t)pos_y;

	return (val_ll*(1.0 - prop_x) + prop_x*val_lr)*(1.0-prop_y) + prop_y*( val_ul*(1.0 - prop_x) + prop_x*val_ur );
}


/// Computes the left-sided difference quatient in x-dim at [it]
real_t Grid::dx_l(const Iterator & it) const {
	return (Cell(it) - Cell(it.Left()))/_geom->Mesh()[0];
}

/// Computes the right-sided difference quatient in x-dim at [it]
real_t Grid::dx_r(const Iterator & it) const {
	return (Cell(it.Right()) - Cell(it)) / _geom->Mesh()[0];
}

/// Computes the left-sided difference quatient in y-dim at [it]
real_t Grid::dy_l(const Iterator & it) const {
	return (Cell(it) - Cell(it.Down())) / _geom->Mesh()[1];
}

/// Computes the right-sided difference quatient in x-dim at [it]
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
	return ((Cell(it.Right()) - 2.0*Cell(it) + Cell(it.Left())) / (_geom->Mesh()[0]) * _geom->Mesh()[0]) ;
}

/// Computes the central difference quatient of 2nd order in y-dim at [it]
real_t Grid::dyy(const Iterator & it) const {
	return ((Cell(it.Top()) - 2.0*Cell(it) + Cell(it.Down())) / (_geom->Mesh()[1]) * _geom->Mesh()[1]);
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
	// Value of u at the right neighbor of the iterator cell (u_{i+1,j})
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
	// Value of v at the upper neighbor (top) of the iterator cell (v_{i,j+1})
	const real_t val_v_t = Cell(it.Top());
	// Value of v at the lower neighbor (down) of the iterator cell (v_{i,j-1})
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
	real_t max = abs(_data[0]);
	const index_t num_cells = _geom->Size()[0] * _geom->Size()[1];
	for (index_t i = 1; i < num_cells; i++) {
		if (max < abs(_data[i])) max = abs(_data[i]);
	}
	return max;
}

/// Returns a pointer to the raw data
real_t * Grid::Data() {
	return _data;
}
