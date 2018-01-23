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


//------------------------------------------------------------------------------
#ifndef __ITERATOR_HPP
#define __ITERATOR_HPP

//------------------------------------------------------------------------------
#include "typedef.hpp"
#include "geometry.hpp"

//------------------------------------------------------------------------------
/** Iterator base class
*/
class Iterator {
public:

  /// Default Constructor
  Iterator();
  /// Constructs a new Iterator depending on a geometry
  Iterator(const Geometry *geom);
  /// Constructs a new Iterator on a geometry with a defined starting value
  Iterator(const Geometry *geom, const index_t &value);

  /// Returns the current position value
  virtual const index_t &Value() const;
  /// Cast operator to convert Iterators to integers
  virtual operator const index_t &() const;
  /// Returns the position coordinates
  virtual multi_index_t Pos() const;

  /// Sets the iterator to the first element
  virtual void First();
  /// Goes to the next element of the iterator, disables it if position is end
  virtual void Next();

  /// Checks if the iterator still has a valid value
  virtual bool Valid() const;

  /// Returns an Iterator that is located left from this one.
  // if we are at the left boundary, the cell sees itself
  virtual Iterator Left() const;

  /// Returns an Iterator that is located right from this one
  // If we are at the right boundary, the cell sees itself
  virtual Iterator Right() const;

  /// Returns an Iterator that is located above this one
  // If we are at the upper domain boundary, the cell sees itself
  virtual Iterator Top() const;

  /// Returns an Iterator that is located below this one
  // If we are at the lower domain boundary, the cell sees itself
  virtual Iterator Down() const;

protected:
  const Geometry *_geom;
  index_t _value;
  bool _valid;
};

//------------------------------------------------------------------------------
/** Iterator for interior cells
*/
class InteriorIterator : public Iterator {
public:
  /// Construct a new InteriorIterator
  InteriorIterator(const Geometry *geom);

  /// Sets the iterator to the first element
  void First();
  /// Goes to the next element of the iterator, disables it if position is end
  void Next();
};

//------------------------------------------------------------------------------
/** Iterator for all cells linewise backwards
*/
class BackwardsIterator : public Iterator {
public:
  /// Construct a new BackwardsIterator
  BackwardsIterator(const Geometry *geom);

  /// Sets the iterator to the last element of the first line
  void First();
  /// Goes to the next element of the iterator, disables it if position is end
  void Next();
};

//------------------------------------------------------------------------------
/** Iterator for domain boundary cells
*/
class BoundaryIterator : public Iterator {
public:
  /// Constructs a new BoundaryIterator
  BoundaryIterator(const Geometry *geom);

  /// Sets the iterator to the first element
  void First();
  /// Goes to the next element of the iterator, disables it if position is end
  void Next();

private:
  index_t _boundary;
};

/------------------------------------------------------------------------------
/** Iterator for interior cells in multigrid algorithm
*/
class MGInteriorIterator : public Iterator {
public:
  /// Construct a new InteriorIterator
  MGInteriorIterator(const Geometry *geom, const MultiGrid *multigrid, index_t cellsize);
  /// Constructs a new BoundaryIterator to start iterating from chosen value
  MGInteriorIterator(const Geometry *geom, const MultiGrid *multigrid, index_t cellsize, index_t value);

  /// Sets the iterator to the first element
  void First();
  /// Goes to the next element of the iterator, disables it if position is end
  void Next();
  /// Returns an Iterator that is located left from this one.
  // if we are at the left boundary, the cell sees itself
  // if there is more than one cell left, we choose the lower one
  MGInteriorIterator Left() const;

  /// Returns an Iterator that is located right from this one
  // If we are at the right boundary, the cell sees itself
  // if there is more than one cell right, we choose the lower one
  MGInteriorIterator Right() const;

  /// Returns an Iterator that is located above this one
  // If we are at the upper domain boundary, the cell sees itself
  // if there is more than one cell at the top, we choose the left one
  MGInteriorIterator Top() const;

  /// Returns an Iterator that is located below this one
  // If we are at the lower domain boundary, the cell sees itself
  // if there is more than one cell below, we choose the left one
  MGInteriorIterator Down() const;
  
private:
    MultiGrid _multigrid;
    index_t _cellsize;
};

//------------------------------------------------------------------------------
/** Iterator for domain boundary cells in multigrid algorithm
*/
class MGBoundaryIterator : public Iterator {
public:
  /// Constructs a new BoundaryIterator
  MGBoundaryIterator(const Geometry *geom, const MultiGrid *multigrid, index_t cellsize);
  /// Constructs a new BoundaryIterator to start iterating from chosen value
  MGBoundaryIterator(const Geometry *geom, const MultiGrid *multigrid, index_t cellsize, index_t value);

  
  /// Sets the iterator to the first element
  void First();
  /// Goes to the next element of the iterator, disables it if position is end
  void Next();  
  /// Returns an Iterator that is located left from this one.
  // if we are at the left boundary, the cell sees itself
  // if there is more than one cell left, we choose the lower one
  MGBoundaryIterator Left() const;

  /// Returns an Iterator that is located right from this one
  // If we are at the right boundary, the cell sees itself
  // if there is more than one cell right, we choose the lower one
  MGBoundaryIterator Right() const;

  /// Returns an Iterator that is located above this one
  // If we are at the upper domain boundary, the cell sees itself
  // if there is more than one cell at the top, we choose the left one
  MGBoundaryIterator Top() const;

  /// Returns an Iterator that is located below this one
  // If we are at the lower domain boundary, the cell sees itself
  // if there is more than one cell below, we choose the left one
  MGBoundaryIterator Down() const;

private:
  index_t _boundary;
  MultiGrid _multigrid;
  index_t _cellsize;
};

//------------------------------------------------------------------------------
#endif // __ITERATOR_HPP
