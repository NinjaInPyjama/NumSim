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
#include "typedef.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "iterator.hpp"
#include <fstream>
#include "grid.hpp"

using namespace std;

int main(int argc, char **argv) {
	const Geometry* geom = new Geometry();
	Iterator it(geom);
    
    const multi_real_t offset = multi_real_t(0.0, 0.0);
    
    Grid * grd = new Grid(geom,offset);
    grd->Initialize(1);
    cout << grd->Interpolate(multi_real_t(2.5, 1.5)) << endl;
    
	//cout << "Value of iterator: " << it.Value() << endl;
	
	//int in = 0;
	//cin >> in;
    /*
    ifstream fin("default.param");
    real_t a;
    string name;
    string gleich;
    while (fin >> name >> gleich >> a){
        cout << name << gleich << a << endl;
    }
    */
    //const Parameter * param = new Parameter();
    //cout << param->IterMax();
    
    //cout << geom->Size()[1] << endl;
    
    
    
    
	return 0;
}
