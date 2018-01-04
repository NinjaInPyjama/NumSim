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
#include <string>

#include "typedef.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "compute.hpp"
#include "visu.hpp"
#include "vtk.hpp"

using namespace std;

int main(int argc, char **argv) {
  // Create parameter and geometry instances with default values
  Parameter param;
  Geometry geom;

  // Create the fluid solver
  Compute comp(&geom, &param);
  
  #ifdef USE_DEBUG_VISU
  // Create and initialize the visualization
  Renderer visu(geom.Length(), geom.Mesh());
  real_t quotient = geom.Length()[1]/geom.Length()[0];
  
  visu.Init(800, 800*quotient);
#endif // USE_DEBUG_VISU

  // Create a VTK generator
 // VTK vtk(geom.Mesh(), geom.Length());

  int32_t start = 0.0;
  
  // Create a VTK File in the folder VTK (must exist)
  
  //vtk.Init("VTK/field");
  //vtk.AddField("Velocity", comp.GetU(), comp.GetV());
  //vtk.AddScalar("Pressure", comp.GetP());
  //vtk.Finish();

  #ifdef USE_DEBUG_VISU
  const Grid *visugrid;

  visugrid = comp.GetVelocity();
#endif // USE_DEBUG_VISU
  
  // Steps
  int i = 0;
  start = clock();
  
  string str(argv[1]);
  
  ofstream myfile;
  do
  {
    myfile.open("output_" + str + ".txt", std::ios_base::app);
  }
  while(!myfile.is_open());
  myfile << param.Re() << " 0 0 0\n";
  
  while (comp.GetTime() <= param.Tend()) {
      
      #ifdef USE_DEBUG_VISU
    // Render and check if window is closed
    switch (visu.Render(visugrid)) {
    case -1:
      return -1;
    case 0:
      visugrid = comp.GetVelocity();
      break;
    case 1:
      visugrid = comp.GetU();
      break;
    case 2:
      visugrid = comp.GetV();
      break;
    case 3:
      visugrid = comp.GetP();
      break;
    case 4:
      visugrid = comp.GetVorticity();
      break;
    default:
      break;
    };
 #endif // USE_DEBUG_VISU
    
	  i++;
	  //cout << "TimeStep: " << i << ", Simulated Time: " << comp.GetTime() << "s / " << param.Tend() << "s" << endl;
	  comp.TimeStep(false);
        
          myfile << comp.GetTime() << " " << comp.GetU()->Cell(120,5) << " " << comp.GetU()->Cell(64,64) << " " << comp.GetU()->Cell(5,120) <<"\n";
	  
	 // vtk.Init("VTK/field");
	 // vtk.AddField("Velocity", comp.GetU(), comp.GetV());
	 // vtk.AddScalar("Pressure", comp.GetP());
	 // vtk.Finish();
	  
  }
  
    
 
  myfile.close();  
	
  cout << "Simulation " << str << " terminated with run time: " << (float)(clock() - start)/1000000 << " s" << endl;

  
 // vtk.Init("VTK/field");
 // vtk.AddField("Velocity", comp.GetU(), comp.GetV());
 // vtk.AddScalar("Pressure", comp.GetP());
 // vtk.Finish();
  
  return 0;
}
