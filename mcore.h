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

using namespace std;
using namespace std::chrono;

/////////////////////////////////////////////////////////
///////////////////// Input /////////////////////////////
/////////////////////////////////////////////////////////

void class_mcube::init() 
{int i,m;
 time_t start, end1;

	time(&start);

	cout << "Range of sectors:" << endl;
	cout << "sec i  j  k "<<endl;	

	for(i=0;i<set;i++){
					nrange(i);																													//range of sectors
					nrange_db(i);}																											//range of sectors debuging

	if(Ismax!=0){
					set_us();																														//set vector for source current us in the mesh		not finnished for X and Z current
					matrix_snA();}																											//average vector potencial for source current

///////////////////////////////////////////////////////////////

	time(&end1);
	cout << "Time: " << difftime(end1, start) << " sec" << endl; 
	cout << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_mcube::cycle_dt()
{int m, i, secstep=1, sum,psum, tmp=0,its=0;
 double h, ht0,l,max,pmax,maxl,t1,t2=0,t2a,t4,t5,t7,t8,k;
 time_t start1, start2, start3, start4, start5, start8, end1, end2, end3, end4, end5, end8;

	Ba=0;Is=0;
	vec_cpy(pBa,zero);
	vec_cpy(B_a,zero);

	if((Bshape==0)||(Bshape==2)) {its=1;}
		else{its=1.0+(ns/4.0);}	

	for(it=its;it<=step;it++){

			stepall=0;

			set_Ba(it);																																					//set instant applied field B_a, Ba according time line
			if(Ismax!=0) 	set_Is(it);																														//set current source Is according time line
			time_variable();																																		//move Jo+=dJ, dJ=0, Ao+=dA, dA=0, T0+=dT, dT=0, set Is,dI,dAs,dA according time line
			set_Jc(it);																																					//calculates Jc according time line, set htpx,htnx/y/z ....

			////////////			chose h 			////////////
			if((htpz<=htpx)&&(htpz<=htpy)&&(ncz!=1)) {h=htpz;}
			if(ncz==1){h=htpz;htpx=0;htpy=0;htnx=0;htny=0;}

			if(ncz!=1){
			if((htpx<=htpy)&&(htpx<=htpz)) {h=htpx;} 
			if((htpy<=htpz)&&(htpy<=htpx)) {h=htpy;}} 

			if(ncx==1){h=htpx;}

			////////////////////////////////////////////////

			ht0=h;

			time(&start3);
			do{
					secstep=1;																													
					max=1e10;
					l=0.9;
					pmax=max;
					maxl=max;

					auto t_start5 = high_resolution_clock::now();

					do{
							if(max<=pmax){pmax = max;}	
							psum=sum;		
							save_temp();																																	//save dJ, dT from previous iteration
							tmp=0;t2=0;t4=0;
							secstep++;

							auto t_start1 = high_resolution_clock::now();

							for(i=0;i<set;i++){							
											recal_U();
											nrange(i);

											#pragma omp parallel for private(m) num_threads(num_threads)		
											for(m=0;m<xsall[i];m++){sector_iter(m);}

											if(sym==2){symm_dT();}	

											if((cube==1)||((cube==0)&&(rel!=1))){damping(l);}	

											cal_dJfromdT();

											if(sym==2){symm_dJ();}

											auto t_start2 = high_resolution_clock::now();													
											cal_av_dA(i);																												//recalculate average potentital according dJ
											auto t_end2 = high_resolution_clock::now();
											t2+=duration<double, milli>(t_end2-t_start2).count();

											if(rel!=1){
																		auto t_start4 = high_resolution_clock::now();
																		magnetic_field();																				//recalculate magnetic field B = B_J(J) + Ba
																		auto t_end4 = high_resolution_clock::now();
																		t4+=duration<double, milli>(t_end4-t_start4).count();

																		if((rel==2)||(rel==3)){
																									recal_theta();														//recalculates unit vector of magnetic field in the cells																	
																									critical_Jc();}														//recalculates critical current density according magnetic field

																		if(nB==1){recal_nB();}}}																//recalculates n power law exponent according magnetic field									

							auto t_end1 = high_resolution_clock::now();
							t1=duration<double, milli>(t_end1-t_start1).count();

							dif(max,1.0,sum);
							if(psum==sum){l=l-0.01;}

							if((cube==1)||((cube==0)&&(rel!=1))){
										maxl = max;
										if(max>pmax){
														do{
																	dif(maxl,l,sum);																					//find maximum difference included dampig factor between two iteration steps 
																	if(maxl>pmax) {l=l-0.01;}} 
															while((maxl>pmax)&&(l>0.02));

															damping(l);																										//apply damping factor to the T vector																							
															tmp=1;

															cal_dJfromdT();
															if(sym==2){symm_dJ();}
															cal_av_dA(2);	

															if(rel!=1){
																					magnetic_field();																	//recalculate magnetic field	
																					if((rel==2)||(rel==3)){
																					recal_theta();																		//recalculates unit vector of magnetic field in the cells																	
																					critical_Jc();}																		//recalculates critical current density according magnetic field
																					if(nB==1){recal_nB();}}}}													//recalculates n power law exponent according magnetic field

							cout<<" max "<<max<<" maxl "<<maxl<<" pmax "<<pmax<<"   changed "<<sum<<" h "<<h<<" l "<<l<<" "<<tmp<<"   Minimi: "<<t1/1000<<"s cal_dA: "<<t2<<"ms B "<< t4 << " ms" << endl;
					}				
//					while(secstep<1);
					while(max>h + h/50.0);																						//final

					auto t_end5 = high_resolution_clock::now();
					t5=duration<double, milli>(t_end5-t_start5).count();

					stepall++;
					cout << stepall << " " << secstep-1 << " htpz " << htpz << " htpx " << htpx << " htpy " << htpy << " h " << h << " ht0*tolJ " << ht0*tolJ << " ht0*tolJ- " << ht0*tolJ+tolJ*ht0/10.0 
							 << "    time " << t5/1000 << "s" <<endl;

					k=10;

					htpz = htpz/k;
					htnz = htnz/k;	
					htpx = htpx/k;
					htnx = htnx/k;	
					htpy = htpy/k;
					htny = htny/k;

					if((htpz<=htpx)&&(htpz<=htpy)&&(ncz!=1)) {h=htpz;}

					if(ncz!=1){
					if((htpx<=htpy)&&(htpx<=htpz)) {h=htpx;} 
					if((htpy<=htpz)&&(htpy<=htpx)) {h=htpy;}} 

					if(ncx==1){h=htpx;}
					if(ncy==1){h=htpy;}
					if(ncz==1){h=htpz;}
			}
//			while(stepall<6);																											//5 = tolJ 1e-5																				
			while(h>tolJ*ht0+tolJ*ht0/10.0);

			time(&end3);

			time(&end1);
			cout << "Time: " << difftime(end3, start3) << " sec " << difftime(end3, start3)/3600 << " h " << endl<<endl;  

			cout << "Saving results... "<< endl;

			time(&start3);

			if((rel!=2)||(rel!=3))magnetic_field();																	//recalculate magnetic field
			if(Btrape==1)external_magnetic_field();																	//recalculate external magnetic field

			interpolation_cell();																										//calculates interpolated J, A, T, m, loss, E at the cells
			recal_U();																															//recalculate U

			save_surface();																													//save variables to surfaces rc J  dJ 

			save_3d();																															//save variables from cells n rc 0 J |J| Jcx A_a B_J B T m
			if(Btrape==1)save_plane();																							//save varaibles from external plane

			loss();																																	//loss per cycel, magnetization M, magnetic moment m  

		  if (ncz==1){save_Tz();}																									//save Tz in form for current lines
					else{
							cal_AV();																												//calculates average of Jx,Jy,Jz,Bx,By,Bz,Tx,Ty,Tz over thickness	
							save_AV();}																											//save average of Jx,Jy,Jz,Bx,By,Bz,Tx,Ty,Tz and Tz for current lines over thickness	

			time(&end3);
			cout << "Time: " << difftime(end3, start3) << " sec" << endl; 
			cout << endl;

	//  J_flux();
	}

	cout << "Total Loss: " << Q << " J"<<endl;
	cout << "Linear Loss: " << Ql << " J"<<endl;
	cout << "Superconducting Loss: " << Qs << " J"<<endl;
	cout << "Total loss from first half of magnetization loop: " << Qh1 << " J"<<endl;
	cout << "Total Loss from entire magnetization loop: " << Qh << " J"<<endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

