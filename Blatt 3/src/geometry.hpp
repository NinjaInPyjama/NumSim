/*
 * Copyright (C) 2015   Malte Brunn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __GEOMETRY_HPP
#define __GEOMETRY_HPP
//------------------------------------------------------------------------------

#define _CRT_SECURE_NO_DEPRECATE
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "iterator.hpp"
#include "communicator.hpp"

//------------------------------------------------------------------------------
class Geometry {
public:
    enum {
    boundaryBottom = 14,
    boundaryLeft = 13,
    boundaryTop = 11,
    boundaryRight = 7,
    cornerTopRight = 3,
    cornerBottomRight = 6,
    cornerBottomLeft = 12,
    cornerTopRight = 9,
    inner = 15
  };
    
  /// Constructs a default geometry:
  // driven cavity with 128 x 128 grid, no-slip boundary conditions
  // as shown below
  //
  //       u=1, v=0
  //    -------------
  //    |           |
  // u=0|           |u=0
  // v=0|           |v=0
  //    |           |
  //    |           |
  //    -------------
  //       u=0, v=0
  Geometry();
  Geometry(const Communicator *comm);

  /// Loads a geometry from a file
  void Load(const char *file);

  void InitializeFlags(Grid * flag, Grid * type, Grid * value) const;

  /// Returns the number of cells in each dimension
  const multi_index_t &Size() const;
  /// Returns the total number of cells in each dimension
  const multi_index_t &TotalSize() const;
  /// Returns the length of the domain
  const multi_real_t &Length() const;
  /// Returns the total length of the domain
  const multi_real_t &TotalLength() const;
  /// Returns the meshwidth
  const multi_real_t &Mesh() const;
  /// Returns the initial velocity
  const multi_real_t &Velocity() const;
  /// Returns the initial pressure
  const real_t &Pressure() const;
  
  /// Returns whether the lower left corner is red or black
  const bool & RedBlack() const;
  
  /// Updates the velocity field u
  void Update_U(Grid *u) const;
  /// Updates the velocity field v
  void Update_V(Grid *v) const;
  /// Updates the pressure field p
  void Update_P(Grid *p) const;

private:
  const Communicator *_comm;

  multi_index_t _size;
  multi_index_t _bsize;
  multi_real_t _length;
  multi_real_t _blength;
  multi_real_t _h;
  bool _redblack;

  multi_real_t _velocity;
  real_t _pressure;

  char * _field;
  
};
//------------------------------------------------------------------------------
#endif // __GEOMETRY_HPP
