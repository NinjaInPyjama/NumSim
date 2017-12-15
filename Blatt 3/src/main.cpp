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
#include "typedef.hpp"
#include "communicator.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "visu.hpp"
#include "vtk.hpp"
//#include "zeitgeist.hpp"

#include <iostream>
#include <sys/stat.h>

int main(int argc, char **argv) {
    
  std::cout << "Hallo1" << std::endl;

  // Create parameter and geometry instances with default values
  Communicator comm(&argc, &argv);
    
  std::cout << "Hallo2" << std::endl;
  //Zeitgeist zeit;
  //zeit.Tic();
  
  
  Parameter param;
  
  std::cout << "Hallo3" << std::endl;
  Geometry geom(&comm);
  std::cout << "Hallo4" << std::endl;
  geom.InitializeFlags();
  
  geom.print('f');
  geom.print('v');
  geom.print('t');
  
  
  // Create the fluid solver
  Compute comp(&geom, &param, &comm);
  
std::cout << "Hallo5" << std::endl;
  
  if (comm.getRank() == 0) {
    // check if folder "VTK" exists
    struct stat info;

    if (stat("VTK", &info) != 0) {
      system("mkdir VTK");
    }
  }
std::cout << "Hallo5.11" << std::endl;
// Create and initialize the visualization
#ifdef USE_DEBUG_VISU
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(800,  800, 
            comm.ThreadDim(), comm.ThreadIdx(), comm.getRank() + 1);
#endif // USE_DEBUG_VISU
std::cout << "Hallo5.12" << std::endl;
  // Create a VTK generator;
  // use offset as the domain shift
  multi_real_t offset;
  offset[0] = comm.ThreadIdx()[0] * (geom.Mesh()[0] * (double)(geom.Size()[0] - 2));
  offset[1] = comm.ThreadIdx()[1] * (geom.Mesh()[1] * (double)(geom.Size()[1] - 2));
  VTK vtk(geom.Mesh(), geom.Length(), geom.TotalLength(), offset, comm.getRank(),
          comm.getSize(), comm.ThreadDim());
std::cout << "Hallo5.13" << std::endl;
#ifdef USE_DEBUG_VISU
  const Grid *visugrid;

  visugrid = comp.GetVelocity();
#endif // USE_DEBUG_VISU
  std::cout << "Hallo5.1" << std::endl;
  // Run the time steps until the end is reached
  //while (comp.GetTime() < param.Tend()) {
      
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

    // Create VTK Files in the folder VTK
    // Note that when using VTK module as it is you first have to write cell
    // information, then call SwitchToPointData(), and then write point data.
//     vtk.Init("VTK/field");
//     vtk.AddRank();
//     vtk.AddCellField("Cell Velocity", comp.GetU(), comp.GetV());
//     vtk.SwitchToPointData();
//     vtk.AddPointField("Velocity", comp.GetU(), comp.GetV());
//     vtk.AddPointScalar("Stream", comp.GetStream());
//     
//     vtk.AddPointScalar("Pressure", comp.GetP());
//     
//     vtk.Finish();

 /* std::cout << "Hallo6" << std::endl;
    // Run a few steps
     for (uint32_t i = 0; i < 9; ++i) {
      */   
 
       comp.TimeStep(true);
       comp.TimeStep(true);
       
//      }
//     bool printOnlyOnMaster = !comm.getRank();
//     comp.TimeStep(printOnlyOnMaster);
//   }
    //zeit.Tac();
    
    //if(comm.getRank()==0) std::cout << "Gesamtzeit: " << zeit.Toc() << std::endl;
  return 0;
}
