/*
 * Copyright (c) 2019 Institute of Electrical Engineering, Slovak Academy of Sciences (authors: Milan Kapolka, Enric Pardo)
 *
 * This file is part of MEMEP3Dtool.
 *
 * MEMEP3Dtool is free software: you can redistribute it and/or modify
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version <http://www.gnu.org/licenses/>.
 *
 * MEMEP3Dtool is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * Details of the modeling tool and the physical and mathematical background can be found in (https://doi.org/10.1016/j.jcp.2017.05.001 ; https://doi.org/10.1088/1361-6668/ab016a ; and https://arxiv.org/abs/1605.09610)
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>
#include <string>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <omp.h>

#include "vector.h"
#include "scalar.h"
#include "declaration.h"
#include "core.h"
#include "mcore.h"


using namespace std;
using namespace std::chrono;

time_t start, middle, end1;

int main(void)
{
time(&start);
cout << endl;
cout << "Copyright (c) 2019 Institute of Electrical Engineering, Slovak Academy of Sciences" << endl;
cout << endl;
cout << "Authors of the code: M. Kapolka and E. Pardo" << endl; 
cout << "Address: Institute of Electrical Engineering, Dubravska 9, 84104 Bratislava, Slovakia" << endl;
cout << "Email address: milan.kapolka@savba.sk and enric.pardo@savba.sk" << endl;
cout << "Version of the code: 0.03" << endl;
cout << "Date: 4.12.2019" << endl;
cout << "WEbpage link to code: https://github.com/epardov/MEMEP3Dtool" << endl;
cout << "GNU General Public License version 3 (GNU GPLv3)" << endl;
cout << endl << "==============================================================" << endl;

cout << endl;
cout << "Time start: " << asctime(localtime(&start)) << endl;

class_mcube mcube; 																						// Declaration	

mcube.clear_files();																					//erase all txt files

mcube.input();																								//load input parameters	
mcube.save_par();																							//save parameters for graphs
mcube.save_nodes();																						//save geometry of nodes
mcube.save_edges();																						//save geometry of edges
mcube.save_cells();																						//save geometry of cells		(set 0 all variable Qo, dQ, dQx...)
mcube.save_surfaces();																				//save geometry of surfaces (set 0 all variable Jo, dJ, dJx...)
mcube.save_matrix();																					//calculates all matrix (q, A, E, H)

time(&middle);
cout << "Time: " << difftime(middle, start) << " sec" << endl; 
cout << endl;

mcube.init();																									//allocated memory for all vcube, minimizing class with coarse mesh
mcube.cycle_dt();																							//minimization	 
time(&end1);
cout << "Time: " << difftime(end1, start) << " sec("<< difftime(end1, start)/3600<<" h)"<< endl;
}

