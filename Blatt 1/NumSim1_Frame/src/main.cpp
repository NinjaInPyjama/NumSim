/*
 *  Copyright (C) 2015   Malte Brunn
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <time.h>

#include "typedef.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "compute.hpp"
#include "vtk.hpp"

using namespace std;

int main(int argc, char **argv) {
  // Create parameter and geometry instances with default values
  Parameter param;
  Geometry geom;

  // Create the fluid solver
  Compute comp(&geom, &param);

  // Create a VTK generator
  VTK vtk(geom.Mesh(), geom.Length());

  int32_t start = 0.0;
  
  // Create a VTK File in the folder VTK (must exist)
  
  vtk.Init("VTK/field");
  vtk.AddField("Velocity", comp.GetU(), comp.GetV());
  vtk.AddScalar("Pressure", comp.GetP());
  vtk.Finish();

  
  // Steps
  int i = 0;
  start = clock();
  while (comp.GetTime() <= param.Tend()) {
	  i++;
	  //cout << "TimeStep: " << i << ", Simulated Time: " << comp.GetTime() << "s / " << param.Tend() << "s" << endl;
	  comp.TimeStep(false);

	  /*
	  vtk.Init("VTK/field");
	  vtk.AddField("Velocity", comp.GetU(), comp.GetV());
	  vtk.AddScalar("Pressure", comp.GetP());
	  vtk.Finish();
	  */
  }
  
  cout << "Run Time: " << (float)(clock() - start)/1000.0 << " s" << endl;

  
  vtk.Init("VTK/field");
  vtk.AddField("Velocity", comp.GetU(), comp.GetV());
  vtk.AddScalar("Pressure", comp.GetP());
  vtk.Finish();
  

  int a = 0;
  cin >> a;
  return 0;
}
