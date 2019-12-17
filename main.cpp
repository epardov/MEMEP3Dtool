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
//mcube.load();																								//load file
}

