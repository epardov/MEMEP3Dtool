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

//set input variables 

void class_cube::input() 
{ double cyl=0,R;

char stmp[50];	
int i,j,k,s;
ifstream inp;
inp.open("input.txt");
  
inp >> stmp >> x;
inp >> stmp >> xl;
inp >> stmp >> y;
inp >> stmp >> z;	
inp >> stmp >> full_matrix; 
inp >> stmp >> nsucx;
inp >> stmp >> nncx;
inp >> stmp >> ncy;
inp >> stmp >> ncz;
inp >> stmp >> n_tapes;
inp >> stmp >> nc_tape;
inp >> stmp >> nc_gap;
inp >> stmp >> d_tape;
inp >> stmp >> d_gap;
inp >> stmp >> elc;
inp >> stmp >> tol_elc;
inp >> stmp >> Bamax;
inp >> stmp >> Bamax1;
inp >> stmp >> Bshape;
inp >> stmp >> Btrape;
inp >> stmp >> Ismax;
inp >> stmp >> rcx_plane;
inp >> stmp >> rcy_plane;
inp >> stmp >> rcz_plane;
inp >> stmp >> x_plane;
inp >> stmp >> y_plane;
inp >> stmp >> z_plane;
inp >> stmp >> ncx_plane;
inp >> stmp >> ncy_plane;
inp >> stmp >> ncz_plane;
inp >> stmp >> theta;
inp >> stmp >> fi;
inp >> stmp >> fi1;
inp >> stmp >> uni;
inp >> stmp >> rel;
inp >> stmp >> nB;
inp >> stmp >> measured_points;
inp >> stmp >> measured_fields;
inp >> stmp >> sym;
inp >> stmp >> Ec;			
inp >> stmp >> Jo;
inp >> stmp >> Jol;
inp >> stmp >> rhoR;
inp >> stmp >> dl;
inp >> stmp >> Jcpa;
inp >> stmp >> Jcpe;
inp >> stmp >> Bo;			
inp >> stmp >> N;	
inp >> stmp >> Nl;
inp >> stmp >> m;
inp >> stmp >> f;		
inp >> stmp >> f1;		
inp >> stmp >> ns;
inp >> stmp >> step;
inp >> stmp >> tolJ;
inp >> stmp >> shape;
inp >> stmp >> num_threads;	

inp.close();

	//constants

	q  =1.60217657e-19;	   																																						//charge of electron
	pi =4.0*atan(1.0);																																								//pi = 3.1415...
	ep0=8.854187817620e-12;      																																			//epsilon
	mi0=4.0*pi*1e-7;						 																																			//mi0		
	c  =(1.0/sqrt(ep0*mi0));       																																		//speed of light

	ncx = nsucx+nncx;

	if ((ncx!=1)&&(ncy!=1)&&(ncz!=1)) {cube=1;}
		else {cube=0;}

	lx=ncx;ly=ncy;lz=ncz;

	cx = x/lx; 
	cy = y/ly; 
	cz = z/lz;

	nscx=9;shiftx1=3;shiftx2=6; if(ncx==1) {nscx=1; shiftx1=0; shiftx2=0;}
	nscy=9;shifty1=3;shifty2=6; if(ncy==1) {nscy=1; shifty1=0; shifty2=0;}
	nscz=9;shiftz1=3;shiftz2=6; if(ncz==1) {nscz=1; shiftz1=0; shiftz2=0;}

  if(shape==4){nscz=ncz;shiftz1=0;shiftz2=0;}


	cal_sectors();																																										//set optimum cells per sector (nscx...) and calculates number of sectors (nsx....)
									 
	set=3;if(sym==0){set=1;}

	nc = ncx*ncy*ncz; 
	nall=pow(nc, 2.0);

	nnx = ncx+1;
	nny = ncy+1;
	nnz = ncz+1;
	nn  = nnx*nny*nnz;

	nsurx = nnx*ncy*ncz; 	
	nsury = ncx*nny*ncz;	
	nsurz = ncx*ncy*nnz;	

	nedgex = ncx*nny*nnz;     
	nedgey = nnx*ncy*nnz;    
	nedgez = nnx*nny*ncz;    

	nsxall  = pow(nsurx,2.0);
	nsyall  = pow(nsury,2.0);
	nsxyall = nsurx * nsury;
	nszall  = pow(nsurz,2.0);

	nxall = nc * nsurx;
	nyall = nc * nsury;
	nzall = nc * nsurz;

	cout << "ncx: " << ncx << endl;
	cout << "ncy: " << ncy << endl;
	cout << "ncz: " << ncz << endl;
	cout << "cx " << cx*1000 << " mm" << endl;
	cout << "cy " << cy*1000 << " mm" << endl;
	cout << "cz " << cz*1000 << " mm" << endl; 
	cout << "nscx " << nscx << endl;
	cout << "nscy " << nscy << endl;
	cout << "nscz " << nscz << endl;
	cout << "shiftx1 " << shiftx1 << endl;
	cout << "shiftx2 " << shiftx2 << endl;
	cout << "shifty1 " << shifty1 << endl;
	cout << "shifty2 " << shifty2 << endl;
	cout << "shiftz1 " << shiftz1 << endl;
	cout << "shiftz1 " << shiftz2 << endl;
	cout << "nsx1 " << nsx[0] << endl;
	cout << "nsy1 " << nsy[0] << endl;
	cout << "nsz1 " << nsz[0] << endl;	
	cout << "nsx2 " << nsx[1] << endl;
	cout << "nsy2 " << nsy[1] << endl;
	cout << "nsz2 " << nsz[1] << endl;
	cout << "nsx3 " << nsx[2] << endl;
	cout << "nsy3 " << nsy[2] << endl;
	cout << "nsz3 " << nsz[2] << endl;	
	cout << "nsall1 " << xsall[0] << endl; 
	cout << "nsall2 " << xsall[1] << endl; 
	cout << "nsall3 " << xsall[2] << endl; 
	cout << "nsall " << nsall << endl;

	cout << "nsurx " << nsurx << endl;
	cout << "nsury " << nsury << endl;
	cout << "nsurz " << nsurz << endl;
	cout << "cube " << cube << endl;
	cout << "threads "<< num_threads << endl;
	cout << "shape " << shape << endl;
	if(shape==4)cout <<"stack of tapes" <<endl;		

	//allocation of memory
	vnode = (struct_node*)malloc(sizeof(struct_node)*nn);

	veX = (struct_edgeX*)malloc(sizeof(struct_edgeX)*(nedgex));
	veY = (struct_edgeY*)malloc(sizeof(struct_edgeY)*(nedgey));
	veZ = (struct_edgeZ*)malloc(sizeof(struct_edgeZ)*(nedgez));

	vsX = (struct_surfaceX*)malloc(sizeof(struct_surfaceX)*(nsurx));
	vsY = (struct_surfaceY*)malloc(sizeof(struct_surfaceY)*(nsury));
	vsZ = (struct_surfaceZ*)malloc(sizeof(struct_surfaceZ)*(nsurz));

	vc = (struct_cell*)malloc(sizeof(struct_cell)*nc);

	if(full_matrix==0){
			cmAx	= (double*)malloc(sizeof(double)*nsurx);
			cmAy	= (double*)malloc(sizeof(double)*nsury);
			cmAz	= (double*)malloc(sizeof(double)*nsurz);

			cmHxyz  = (double*)malloc(sizeof(double)*nsury);
			cmHxzy  = (double*)malloc(sizeof(double)*nsurz);

			cmHyzx  = (double*)malloc(sizeof(double)*nsurz);
			cmHyxz  = (double*)malloc(sizeof(double)*nsurx);

			cmHzxy  = (double*)malloc(sizeof(double)*nsurx);
			cmHzyx  = (double*)malloc(sizeof(double)*nsury);}

	if(full_matrix==1){
			mAx = (double**)malloc(sizeof(double*)*nsurx);
			mAy = (double**)malloc(sizeof(double*)*nsury);
			mAz = (double**)malloc(sizeof(double*)*nsurz);

			mHxyz  = (double**)malloc(sizeof(double*)*nc);
			mHxzy  = (double**)malloc(sizeof(double*)*nc);

			mHyzx  = (double**)malloc(sizeof(double*)*nc);
			mHyxz  = (double**)malloc(sizeof(double*)*nc);

			mHzxy  = (double**)malloc(sizeof(double*)*nc);
			mHzyx  = (double**)malloc(sizeof(double*)*nc);

			for(i=0;i<nsurx;i++)  {mAx[i]  = (double*)malloc(sizeof(double)*(nsurx-i));}
			for(i=0;i<nsury;i++)  {mAy[i]  = (double*)malloc(sizeof(double)*(nsury-i));}
			for(i=0;i<nsurz;i++)	{mAz[i]  = (double*)malloc(sizeof(double)*(nsurz-i));}

			for(i=0;i<nc;i++)     {mHyxz[i] = (double*)malloc(sizeof(double)*nsurx);
														 mHzxy[i] = (double*)malloc(sizeof(double)*nsurx);
														 mHxyz[i] = (double*)malloc(sizeof(double)*nsury);
														 mHzyx[i] = (double*)malloc(sizeof(double)*nsury);
														 mHxzy[i] = (double*)malloc(sizeof(double)*nsurz);
														 mHyzx[i] = (double*)malloc(sizeof(double)*nsurz);}}

	if(Ismax!=0){
			snAx	= (double*)malloc(sizeof(double)*nsurx);
			snAy	= (double*)malloc(sizeof(double)*nsury);
			snAz	= (double*)malloc(sizeof(double)*nsurz);}

	Jx	= (double*)malloc(sizeof(double)*(ncx*ncy));
	Jy	= (double*)malloc(sizeof(double)*(ncx*ncy));
	Jz	= (double*)malloc(sizeof(double)*(ncx*ncy));

	Bx	= (double*)malloc(sizeof(double)*(ncx*ncy));
	By	= (double*)malloc(sizeof(double)*(ncx*ncy));
	Bz	= (double*)malloc(sizeof(double)*(ncx*ncy));

	Tx	= (double*)malloc(sizeof(double)*(ncx*ncy));
	Ty	= (double*)malloc(sizeof(double)*(ncx*ncy));
	Tz	= (double*)malloc(sizeof(double)*(ncx*ncy));

	cii	= (int*)malloc(sizeof(int)*nsall);
	cai	= (int*)malloc(sizeof(int)*nsall);
	cij	= (int*)malloc(sizeof(int)*nsall);
	caj	= (int*)malloc(sizeof(int)*nsall);
	cik	= (int*)malloc(sizeof(int)*nsall);
	cak	= (int*)malloc(sizeof(int)*nsall);

	if(Btrape==1){
			nc_plane = ncx_plane * ncy_plane * ncz_plane;
			plx=ncx_plane;
			ply=ncy_plane;
			plz=ncz_plane;

			pcx[0] = x_plane/plx; 
			pcx[1] = y_plane/ply; 
			pcx[2] = z_plane/plz;

			V_ave_plane = pcx[0]*pcx[1]*pcx[2];

			RC_plane = (double**)malloc(sizeof(double*)*nc_plane);
			BJ_plane = (double**)malloc(sizeof(double*)*nc_plane);
			B_plane	 = (double**)malloc(sizeof(double*)*nc_plane);

			emHxyz = (double**)malloc(sizeof(double*)*nc_plane);
			emHxzy = (double**)malloc(sizeof(double*)*nc_plane);
			emHyzx = (double**)malloc(sizeof(double*)*nc_plane);
			emHyxz = (double**)malloc(sizeof(double*)*nc_plane);
			emHzxy = (double**)malloc(sizeof(double*)*nc_plane);
			emHzyx = (double**)malloc(sizeof(double*)*nc_plane);

			for(i=0;i<nc_plane;i++){		
						RC_plane[i] = (double*)malloc(sizeof(double)*3);
						BJ_plane[i] = (double*)malloc(sizeof(double)*3);
						B_plane[i]  = (double*)malloc(sizeof(double)*3);

						emHyxz[i] = (double*)malloc(sizeof(double)*nsurx);
						emHzxy[i] = (double*)malloc(sizeof(double)*nsurx);

						emHxyz[i] = (double*)malloc(sizeof(double)*nsury);
						emHzyx[i] = (double*)malloc(sizeof(double)*nsury);

						emHxzy[i] = (double*)malloc(sizeof(double)*nsurz);
						emHyzx[i] = (double*)malloc(sizeof(double)*nsurz);}}

	if(rel==3){	data = (double***)malloc(sizeof(double**)*measured_fields);

							for(i=0;i<measured_fields;i++){data[i] = (double**)malloc(sizeof(double*)*measured_points);}

							for(i=0;i<measured_fields;i++){
									for(j=0;j<measured_points;j++) {data[i][j]	= (double*)malloc(sizeof(double)*4);}}}//[B][angle][0-point, 1-B, 2-Jc, 3-n]

	if(rel==3){	data1 = (double**)malloc(sizeof(double*)*measured_points);
							for(i=0;i<measured_points;i++) {data1[i] = (double*)malloc(sizeof(double)*4);}} 			//[point][0-point, 1-B, 2-Jc, 3-n]

	//constants

	T=1.0/f;										 																																			//period [s]	
	dt=T/ns;																																													//delta t [s]

	V_ave=x*y*z/(ncx*ncy*ncz);	 																																			//average volume
	V_tot=x*y*z;								 																																			//total volume		
	//Uo=Ec*sqrt(Jcpe*Jcpe+Jcpa*Jcpa)/(N+1.0);						 																						//Uo parameter for anisotropic equation
	Uo=Ec*Jcpe/(N+1.0);						 																																		//Uo parameter for anisotropic equation
	mo=(N+1.0)/2.0;							 																																			//mo factor finterJ_polar(or anisotropic equation	 

	tola=1e-16;									 																																			//tolerance of walls

	Px=0;																																															//loss in dt
	Pix=0;																																														//loss in previous dt
	Q=0;																																															//total loss of sample within cycle
	Qh=0;												 																																			//total loss calculated from hysterezis loop	
	Qh1=0;
	Qs=0;
	Ql=0;
	Plin=0;
	Pilin=0;
	Psup=0;
	Pisup=0;
	rho=0;
	zero[0]=0;
	zero[1]=0;
	zero[2]=0;

	if(nncx!=0){
				rho=dl*rhoR/(nncx*cx); 																																			//real resistivity of linear material
				cout << "rho: " << rho  << endl;}

	if((Jol==0)&&(rho!=0)){
				Jol=(Ec/rho);}         																																			//current density of linear material

	cout << "Jo: "  << Jo	 << endl;
	cout << "Jol: " << Jol << endl;

	Ba = Bamax * sin (2.0 * pi * f * T/4.0);

	B_a[0] = Ba * sin(fi*pi/180.0)*cos(theta*pi/180.0);			if(fabs(B_a[0])<1e-12) B_a[0]=0;  				//instant magnetic field B_a x component		
	B_a[1] = Ba * sin(fi*pi/180.0)*sin(theta*pi/180.0);			if(fabs(B_a[1])<1e-12) B_a[1]=0; 					//instant magnetic field B_a y component
	B_a[2] = Ba * cos(fi*pi/180.0);													if(fabs(B_a[2])<1e-12) B_a[2]=0;  				//instant magnetic field B_a z component

	vec_divi(B_a, Ba, ea); 																																						//eq unit vector of applied field at the peak ea=B_a/Ba;				

	Ba=0;
	vec_cpy(B_a,zero);
	d_ncx=ncx,d_ncy=ncy;


	if(ncz==1){
			ratio=(cy/cx)/(d_ncy/d_ncx);
			cout << "Cartezian ratio " << ratio << endl;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// Geometry nodes ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::node_rc(int i, int j, int k, double XN[6])
{
	
	if(uni==1){																																		//uniform mesh in cartezian coordinate system
				XN[0] = cx * i;
				XN[1] = cy * j;
				XN[2] = cz * k;

				if(shape==4) {XN[2] = k * (d_tape/nc_tape) - k/(nc_tape+1)*(d_tape/nc_tape) + k/(nc_tape+1)*d_gap;}}

	if(uni==2){																																		//semi-uniform mesh
				if ((x==0)||(i==0)) {XN[0]=0;}
				else{						
							 if (i<=nscx/2){
										 XN[0] = ((x-xl)/nscx)*i;}

							 if ((i>nscx/2)&&(i<=(nscx/2) + nncx)){
										 XN[0] = ((x-xl)/2) + (i-(nscx/2))*(xl/nncx);}

						   if (i>(nscx/2) + nncx){
										 XN[0] = ((x-xl)/2) + xl + (i - (nscx/2 + nncx))*((x-xl)/nscx);
										 if (i==ncx) {XN[0]=x;}}}

				if ((y==0)||(j==0))	{XN[1]=0;}else{XN[1]=(y/ncy)*j;}

				if ((z==0)||(k==0))	{XN[2]=0;}else{XN[2]=(z/ncz)*k;}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int class_cube::inode(int i,int j,int k)
{return i + nnx*j + nnx*nny*k;}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::nnodetoijk(int n, int& i, int& j, int& k)
{
	// n = i + nnx*j + nnx*nny*k;
	k=  n/(nnx*nny);
	j= (n - k*nnx*nny)/nnx ;
	i= (n - k*nnx*nny) - (j*nnx);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////// Geometry edges //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::edgeX_rc(int i, int j, int k, double RC[6]) 
{int l=0;

	for(l=0;l<=2;l++)		{RC[l] = (xnode_rc(i, j, k, l) + xnode_rc(i+1, j, k, l)) * 0.5;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::edgeY_rc(int i, int j, int k, double RC[6]) 
{	int l=0;

	for(l=0;l<=2;l++) 	{RC[l] = (xnode_rc(i, j, k, l) + xnode_rc(i, j+1, k, l)) * 0.5;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::edgeZ_rc(int i, int j, int k, double RC[6]) 
{int l=0;

	for(l=0;l<=2;l++)		{RC[l] = (xnode_rc(i, j, k, l) + xnode_rc(i, j, k+1, l)) * 0.5;}									
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::edgeX_lenght(int i, int j, int k)
{return xnode_rc(i+1, j, k, 0) - xnode_rc(i, j, k, 0);}																																	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::edgeY_lenght(int i, int j, int k)
{return xnode_rc(i, j+1, k, 1) - xnode_rc(i, j, k, 1);}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
							 
double class_cube::edgeZ_lenght(int i, int j, int k)
{return xnode_rc(i, j, k+1, 2) - xnode_rc(i, j, k, 2);}																																									

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int class_cube::ijktonedgeX(int i,int j,int k)
{return i + ncx*j + ncx*nny*k;}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int class_cube::ijktonedgeY(int i,int j,int k)
{return i + nnx*j + nnx*ncy*k;}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int class_cube::ijktonedgeZ(int i,int j,int k)
{return i + nnx*j + nnx*nny*k;}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::nedgeXtoijk(int n, int& ix, int& iy, int& iz)
{
	iz= n / (ncx*nny);
	iy= ( n - iz*ncx*nny )/ncx;
	ix= n - iz*ncx*nny - iy*ncx;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::nedgeYtoijk(int n, int& ix, int& iy, int& iz)
{
	iz= n / (nnx*ncy);
	iy= ( n - iz*nnx*ncy )/nnx;
	ix= n - iz*nnx*ncy - iy*nnx;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::nedgeZtoijk(int n, int& ix, int& iy, int& iz)
{
	iz= n / (nnx*nny);
	iy= ( n - iz*nnx*nny )/nnx;
	ix= n - iz*nnx*nny - iy*nnx;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::edgeZtosurXY(int m, int& x1, int& x2, int& y1, int& y2)
{int i=0,j=0,k=0;

  i=veZ[m].adr[0];
	j=veZ[m].adr[1];
	k=veZ[m].adr[2];

	x1 = ijktonsurX(i, j-1, k);
	x2 = ijktonsurX(i, j, k);
	y1 = ijktonsurY(i, j, k);
	y2 = ijktonsurY(i-1, j, k);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::edgeXtosurYZ(int m, int& y1, int& y2, int& z1, int& z2)
{int i=0,j=0,k=0;

  i=veX[m].adr[0];
	j=veX[m].adr[1];
	k=veX[m].adr[2];

	y1 = ijktonsurY(i, j, k-1);
	y2 = ijktonsurY(i, j, k);
	z1 = ijktonsurZ(i, j-1, k);
	z2 = ijktonsurZ(i, j, k);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::edgeYtosurXZ(int m, int& x1, int& x2, int& z1, int& z2)
{int i=0,j=0,k=0;

  i=veY[m].adr[0];
	j=veY[m].adr[1];
	k=veY[m].adr[2];

	x1 = ijktonsurX(i, j, k-1);
	x2 = ijktonsurX(i, j, k);
	z1 = ijktonsurZ(i, j, k);
	z2 = ijktonsurZ(i-1, j, k);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::boundary_addres_edge()
{int i,j,k,m;

	j=0;k=0;

	if(ncz==1){
	for(i=1;i<ncx;i++){////// the first bounday surface //////
				m = ijktonedgeZ(i, j, k);

				veZ[m].adrs[0] = ijktonsurX(i, ncy-1, k);	//address of the surface in the opposite site
				veZ[m].adrs[1] = ijktonsurX(i, j, k); 
				veZ[m].adrs[2] = ijktonsurY(i, j, k);
				veZ[m].adrs[3] = ijktonsurY(i-1, j, k);

				veZ[m].adrc[0] = icell(i, ncy-1, k);
				veZ[m].adrc[1] = icell(i,j,k);
				veZ[m].adrc[2] = icell(i-1, j,k);
				veZ[m].adrc[3] = icell(i-1,ncy-1,k);}


	for(i=1;i<ncx;i++){////// the second bounday surface //////
				m = ijktonedgeZ(i, ncy, 0);

				veZ[m].adrs[0] = ijktonsurX(i, ncy-1, k);	//address of the surface in the opposite site
				veZ[m].adrs[1] = ijktonsurX(i, 0, k); 
				veZ[m].adrs[2] = ijktonsurY(i, ncy, k);
				veZ[m].adrs[3] = ijktonsurY(i-1, ncy, k);

				veZ[m].adrc[0] = icell(i, ncy-1, k);
				veZ[m].adrc[1] = icell(i,0,k);
				veZ[m].adrc[2] = icell(i-1, 0,k);
				veZ[m].adrc[3] = icell(i-1,ncy-1,k);}}

	i=0;j=0;

	if(ncx==1){
	for(k=1;k<ncz;k++){////// the first bounday surface //////
				m = ijktonedgeX(i, j, k);

				veX[m].adrs[0] = ijktonsurY(i, j, k-1);	//address of the surface in the opposite site
				veX[m].adrs[1] = ijktonsurY(i, j, k); 
				veX[m].adrs[2] = ijktonsurZ(i, ncy-1, k);
				veX[m].adrs[3] = ijktonsurZ(i, j, k);

				veX[m].adrc[0] = icell(i, j, k-1);
				veX[m].adrc[1] = icell(i,j,k);
				veX[m].adrc[2] = icell(i, ncy-1,k);
				veX[m].adrc[3] = icell(i,ncy-1,k-1);}

	j=ncy;

	for(k=1;k<ncz;k++){////// the second bounday surface //////
				m = ijktonedgeX(i, j, k);

				veX[m].adrs[0] = ijktonsurY(i, j, k-1);	//address of the surface in the opposite site
				veX[m].adrs[1] = ijktonsurY(i, j, k); 
				veX[m].adrs[2] = ijktonsurZ(i, j-1, k);
				veX[m].adrs[3] = ijktonsurZ(i, 0, k);

				veX[m].adrc[0] = icell(i,0, k-1);
				veX[m].adrc[1] = icell(i,0,k);
				veX[m].adrc[2] = icell(i, j-1,k);
				veX[m].adrc[3] = icell(i,j-1,k-1);}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////// Geometry surfaces /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::surfacex_rc(int i, int j, int k, double XSC[6])	
{int l;

	for(l=0;l<=2;l++)	{XSC[l] = (xnode_rc(i, j, k, l) + xnode_rc(i, j, k+1, l) + xnode_rc(i, j+1, k, l) + xnode_rc(i, j+1, k+1, l)) * 0.25;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::surfacey_rc(int i, int j, int k, double YSC[6])	
{int l;

	for(l=0;l<=2;l++)	{YSC[l] = (xnode_rc(i, j, k, l) + xnode_rc(i, j, k+1, l) + xnode_rc(i+1, j, k, l) + xnode_rc(i+1, j, k+1, l)) * 0.25;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::surfacez_rc(int i, int j, int k, double ZSC[6])	
{int l;

	for(l=0;l<=2;l++)	{ZSC[l] = (xnode_rc(i, j, k, l) + xnode_rc(i+1, j, k, l) + xnode_rc(i, j+1, k, l) + xnode_rc(i+1, j+1, k, l)) * 0.25;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::surfaceX_size(int i, int j, int k, double SS[4])
{

	if (i==0)							{SS[0] = xnode_rc(i+1,j,k,0)  - xnode_rc(i,j,k,0);}
	if ((i!=0)&&(i!=ncx)) {SS[0] = (1.0/2.0)*(xnode_rc(i+1,j,k,0)  - xnode_rc(i,j,k,0)) + (1.0/2.0)*(xnode_rc(i,j,k,0) - xnode_rc(i-1,j,k,0));}				//a
	if (i==ncx)						{SS[0] = xnode_rc(i,j,k,0) 	 - xnode_rc(i-1,j,k,0);}

												SS[1] = xnode_rc(i,j+1,k,1) - xnode_rc(i,j,k,1);																																						//b		
												SS[2] = xnode_rc(i,j,k+1,2) - xnode_rc(i,j,k,2);																																						//c
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			

void class_cube::surfaceY_size(int i, int j, int k, double SS[4])
{	
	SS[0] = xnode_rc(i+1,j,k,0) - xnode_rc(i,j,k,0);																																																	//a
	SS[2] = xnode_rc(i,j,k+1,2) - xnode_rc(i,j,k,2);																																																	//c

	if (j==0)							{SS[1] = xnode_rc(i,j+1,k,1) - xnode_rc(i,j,k,1);}																																					//b
	if ((j!=0)&&(j!=ncy)) {SS[1] = (1.0/2.0)*(xnode_rc(i,j+1,k,1) - xnode_rc(i,j,k,1)) + (1.0/2.0)*(xnode_rc(i,j,k,1) - xnode_rc(i,j-1,k,1));}
	if (j==ncy)           {SS[1] = xnode_rc(i,j,k,1) - xnode_rc(i,j-1,k,1);}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::surfaceZ_size(int i, int j, int k, double SS[4])
{	
	SS[0] = xnode_rc(i+1,j,k,0) - xnode_rc(i,j,k,0);																																																	//a
	SS[1] = xnode_rc(i,j+1,k,1) - xnode_rc(i,j,k,1);																																																	//b

	if (k==0)							{SS[2] = xnode_rc(i,j,k+1,2) - xnode_rc(i,j,k,2);}
	if ((k!=0)&&(k!=ncz)) {SS[2] = (1.0/2.0)*(xnode_rc(i,j,k+1,2) - xnode_rc(i,j,k,2)) + (1.0/2.0)*(xnode_rc(i,j,k,2) - xnode_rc(i,j,k-1,2));}				//c
	if (k==ncz) 					{SS[2] = xnode_rc(i,j,k,2) - xnode_rc(i,j,k-1,2);} 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::surX_S(int m, int i, int j, int k)
{return vsX[m].size[1] * vsX[m].size[2];}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
double class_cube::surY_S(int m, int i, int j, int k)
{return vsY[m].size[0] * vsY[m].size[2];}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::surZ_S(int m, int i, int j, int k)
{return vsZ[m].size[0] * vsZ[m].size[1];}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::vol_infX(int i, int j, int k)
{int m=0;
 double V=0;

	m = icell(i,j,k);	 

	if (i==0) 						{return vc[m].V;} 
	if ((i!=0)&&(i!=ncx)) {return (1.0/2.0)*vc[m].V + (1.0/2.0)*vc[m-1].V;}
	if (i==ncx) 					{return vc[m-1].V;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::vol_infY(int i, int j, int k)
{int m=0;

	m = icell(i,j,k);	

	if (j==0) 						{return vc[m].V;}
	if ((j!=0)&&(j!=ncy)) {return (1.0/2.0)*vc[m].V + (1.0/2.0)*vc[m - ncx].V;}
	if (j==ncy)						{return vc[m - ncx].V;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::vol_infZ(int i, int j, int k)
{int m=0;

	m = icell(i,j,k);	

	if (k==0)							{return vc[m].V;}
	if ((k!=0)&&(k!=ncz)) {return (1.0/2.0)*vc[m].V + (1.0/2.0)*vc[m - ncx*ncy].V;}
	if (k==ncz) 					{return vc[m - ncx*ncy].V;}
}		

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::surf_shape()
{int i=0;
 double S[3]={0}, R=0, W[3]={0};

	if ((shape==1)||(shape==2)){
											S[0]=x/2;
											S[1]=y/2;
											if(shape==1) {S[2]=z/2;}

							//			R=sqrt(pow(S[0],2.0))-(cx/4);
											R=sqrt(pow(S[0],2.0));
											for(i=0;i<nc;i++){
														vec_subs2(vc[i].rc,S,W);
														if(shape==2) W[2]=0;

														if(vec_mod(W)>R){
																vc[i].N  = Nl;
																vc[i].Jc	=	Jol;}}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int class_cube::ijktonsurX(int i,int j,int k)
{return i + nnx*j + nnx*ncy*k;}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int class_cube::ijktonsurY(int i,int j,int k)
{return i + ncx*j + nny*ncx*k;}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int class_cube::ijktonsurZ(int i,int j,int k)
{return i + ncx*j + ncx*ncy*k;}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::nsurXtoijk(int n, int& ix, int& iy, int& iz)
{
	iz= n / (nnx*ncy);
	iy= ( n - iz*nnx*ncy )/nnx;
	ix= ( n - iz*nnx*ncy - iy*nnx );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::nsurYtoijk(int n, int& ix, int& iy, int& iz)
{
	iz= n / (nny*ncx);
	iy= ( n - iz*nny*ncx )/ncx;
	ix= ( n - iz*nny*ncx - iy*ncx);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::nsurZtoijk(int n, int& ix, int& iy, int& iz)
{
	iz= n / (ncy*ncx);
	iy= ( n - iz*ncy*ncx )/ncx;
	ix= ( n - iz*ncy*ncx - iy*ncx );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Geometry cell /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::xnode_rc(int i, int j, int k, int l)
{double r;

	if (uni==1){ 																										//uniform mesh in cartezian coordinate system
							 if(l==0)	{return cx * i;} 											 	
							 if(l==1)	{return cy * j;}
							 if(l==2)	{r = cz * k;
											 	 if(shape==4){r = k * (d_tape/nc_tape) - k/(nc_tape+1)*(d_tape/nc_tape) + k/(nc_tape+1)*d_gap;}
												 return r;}}

	if (uni==2){																																	//linear material non-uniform mesh
								if(l==0){
												if ((x==0)||(i==0)) r=0;
														else	
														 		{						
																	 if (i<=nscx/2)													{r = ((x-xl)/nscx)*i;}
																	 if ((i>nscx/2)&&(i<=(nscx/2) + nncx))	{r = ((x-xl)/2) + (i-(nscx/2))*(xl/nncx);}
																	 if (i>(nscx/2) + nncx)									{r = ((x-xl)/2) + xl + (i - (nscx/2 + nncx))*((x-xl)/nscx);
																				 																	 if (i==ncx) {r=x;}}}
										 		return r;}
								if(l==1){return (y/ncy)*j;}
								if(l==2){return (z/ncz)*k;}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::xnode_prc(int i, int j, int k, int l)
{
		if(l==0){	if(ncx_plane==1) {return x_plane;}
								else{return rcx_plane - x_plane/2.0 + pcx[0] * i;}}					
		if(l==1){	if(ncy_plane==1) {return y_plane;}
								else{return rcy_plane - y_plane/2.0 + pcx[1] * j;}}
		if(l==2){return rcz_plane - pcx[2]/2.0 + pcx[2] * k;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::cell_rc(int i, int j, int k, double XC[6])	
{int l; 

	for(l=0;l<=2;l++){
					XC[l] = (xnode_rc(i, j, k, l)     + xnode_rc(i+1, j, k, l) + xnode_rc(i, j+1, k, l) + xnode_rc(i+1, j+1, k, l) + xnode_rc(i, j, k+1, l) + xnode_rc(i+1, j, k+1, l) 
								 + xnode_rc(i, j+1, k+1, l) + xnode_rc(i+1, j+1, k+1, l)) * 0.125;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::cell_rc_plane(int i, int j, int k, int l)	
{ 
 	return (xnode_prc(i, j, k, l)     + xnode_prc(i+1, j, k, l) + xnode_prc(i, j+1, k, l) 
		    + xnode_prc(i+1, j+1, k, l) + xnode_prc(i, j, k+1, l) + xnode_prc(i+1, j, k+1, l) 
				+ xnode_prc(i, j+1, k+1, l) + xnode_prc(i+1, j+1, k+1, l)) * 0.125;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::cell_vol(double SC[3], int i, int j, int k)  	
{return vol(SC);}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::cell_size(int i, int j, int k, double SC[3])
{ 
	SC[0] = xnode_rc(i+1,j,k,0) - xnode_rc(i,j,k,0);																																					//dx
	SC[1] = xnode_rc(i,j+1,k,1) - xnode_rc(i,j,k,1);																																					//dy
	SC[2] = xnode_rc(i,j,k+1,2) - xnode_rc(i,j,k,2);																																					//dz
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::cell_lin ()
{int i=0, j=0, k=0, n=0,a=0;

	if(shape==0){
			for (j=0;j<ncy;j++){
						for (i=0;i<nncx;i++){
										n=(nsucx/2) + i + j*ncx;
										vc[n].Jc=Jol;		
										vc[n].N = Nl;}}}

	if(shape==3){
			for (j=0;j<ncy;j++){
						for (i=0;i<nncx;i++){
										n=(nsucx/(nncx+1)) + i*(nsucx/(nncx+1))+i + j*ncx;
										vc[n].Jc=Jol;
										vc[n].N = Nl;}}}

	if(shape==4){
			for(i=1;i<n_tapes;i++){
					if(i>1)a=1;else{a=0;}
							for(n=0;n<nc;n++){
									if(vc[n].adr[2]== i*nc_tape + (i-1)*a){
									vc[n].Jc=Jol;
									vc[n].N = Nl;}}}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int class_cube::icell(int i,int j,int k)
{return i + ncx*j + ncx*ncy*k;}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int class_cube::iplanecell(int i,int j,int k)
{
	return i + ncx_plane*j + ncx_plane*ncy_plane*k;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::ncelltoijk(int n, int& i, int& j, int& k)
{
	// n = i + ncx*j + nx*ny*k;
	k=  n/(ncx*ncy);
	j= (n - k*ncx*ncy)/ncx ;
	i= (n - k*ncx*ncy) - (j*ncx);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::ncellplanetoijk(int n, int& i, int& j, int& k)
{
	// n = i + ncx*j + nx*ny*k;
	k=  n/(ncx_plane*ncy_plane);
	j= (n - k*ncx_plane*ncy_plane)/ncx_plane;
	i= (n - k*ncx_plane*ncy_plane) - (j*ncx_plane);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Geometry sector ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::nrange(int set)
{int n,i,j,k,m;

		/////////////// range of set 1 cells /////////////////////////////
		if(set==0){
				for(m=0;m<xsall[set];m++){
								nsectortoijk(set,m,i,j,k);

								if(ncx==1){cii[m]=0;cai[m]=0;}
										else{
													cii[m] = nscx*i;
													if(i==0){cai[m]=nscx-1;}
														else{
																	cai[m]=cii[m]+nscx-1;
																	if((i==nsx[0]-1)&&(sym==1)){cai[m]=NCX-1;}}}

								if(ncy==1){cij[m]=0;caj[m]=0;}
										else{	
													cij[m] = nscy*j;	
													if(j==0){caj[m]=nscy-1;}
														else{
																	caj[m]=cij[m]+nscy-1;
																	if((j==nsy[0]-1)&&(sym==1)){caj[m]=NCY-1;}}}

								if(ncz==1){cik[m]=0;cak[m]=0;}
										else{
													cik[m] = nscz*k;
													if(k==0){cak[m]=nscz-1;}
															else{
																	cak[m]=cik[m]+nscz-1;
																	if((k==nsz[0]-1)&&(sym==1)){cak[m]=NCZ-1;}}}}}

		/////////////// range of set 2 cells /////////////////////////////
		if(set==1){
				for(m=0;m<xsall[set];m++){
								nsectortoijk(set,m,i,j,k);

								if(ncx==1){cii[m]=0;cai[m]=0;}
										else{	
													if(i==0){cii[m] = 0;}else{cii[m] = nscx*i - (nscx-shiftx1);}
													if(i==0){cai[m] = shiftx1-1;}else{cai[m] = cii[m] + nscx - 1;}
													if((i==nsx[1]-1)&&(sym==1)){cai[m] = NCX-1;}}

								if(ncy==1){cij[m]=0;caj[m]=0;}
										else{	
													if(j==0){cij[m] = 0;}else{cij[m] = nscy*j - (nscy-shifty1);}
													if(j==0){caj[m] = shifty1-1;}else{caj[m] = cij[m] + nscy - 1;}
													if((j==nsy[1]-1)&&(sym==1)){caj[m] = NCY-1;}}

								if(ncz==1){cik[m]=0;cak[m]=0;}
										else{		
													if(k==0){cik[m] = 0;}else{cik[m] = nscz*k - (nscz-shiftz1);}
													if(k==0){cak[m] = shiftz1-1;}else{cak[m] = cik[m] + nscz - 1;}
													if((k==nsz[1]-1)&&(sym==1)){cak[m] = NCZ-1;}}}}

		/////////////// range of set 3 cells /////////////////////////////
		if(set==2){
				for(m=0;m<xsall[set];m++){
								nsectortoijk(set,m,i,j,k);
								if(ncx==1){cii[m]=0;cai[m]=0;}
										else{	
													if(i==0){cii[m] = 0;}else{cii[m] = nscx*i - (nscx-shiftx2);}
													if(i==0){cai[m] = shiftx2-1;}else{cai[m] = cii[m] + nscx - 1;}
													if((i==nsx[2]-1)&&(sym==1)){cai[m] = NCX-1;}}

								if(ncy==1){cij[m]=0;caj[m]=0;}
										else{	
													if(j==0){cij[m] = 0;}else{cij[m] = nscy*j - (nscy-shifty2);}
													if(j==0){caj[m] = shifty2-1;}else{caj[m] = cij[m] + nscy - 1;}
													if((j==nsy[2]-1)&&(sym==1)){caj[m] = NCY-1;}}



								if(ncz==1){cik[m]=0;cak[m]=0;}
										else{		
													if(k==0){cik[m] = 0;}else{cik[m] = nscz*k - (nscz-shiftz2);}
													if(k==0){cak[m] = shiftz2-1;}else{cak[m] = cik[m] + nscz - 1;}
													if((k==nsz[2]-1)&&(sym==1)){cak[m] = NCZ-1;}}}}

	if(sym==0){
			cii[0]=0;
			cai[0]=ncx-1;
			cij[0]=0;
			caj[0]=ncy-1;
			cik[0]=0;
			cak[0]=ncz-1;}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::nrange_db(int set)
{int i,j,k,m;
	
	for(m=0;m<xsall[set];m++){
							nsectortoijk(set,m, i, j, k);
							cout << "set " << set << " s " << m << " " << i << " " << j << " " << k << " 	" << cii[m] << " " << cai[m] <<  " " << cij[m] << " " << caj[m] << " " << cik[m] << " " << cak[m] <<endl;}
	cout << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::cal_sectors()
{int lcs, sector;

	///set optimum///
	sector=1;

	//////////////////////// nsx ///////////////////////////
	if(sym==1){NCX=ncx;}
		else{NCX=(ncx/2)+sector;}

	if(ncx==1){nsx[0]=1;nsx[1]=nsx[0];nsx[2]=nsx[0];}
		else{
				nsx[0] = NCX/nscx;					
				lcs    = NCX%nscx;									if(lcs>0){nsx[0]++;}
				nsx[1] = (NCX-shiftx1)/nscx;				nsx[1]++;			
				lcs    = (NCX-shiftx1)%nscx;				if(lcs>0){nsx[1]++;}
				nsx[2] = (NCX-shiftx2)/nscx;				nsx[2]++;			
				lcs    = (NCX-shiftx2)%nscx;				if(lcs>0){nsx[2]++;}}

	//////////////////////// nsy ///////////////////////////
	if(sym==1){NCY=ncy;}
		else{NCY=(ncy/2)+sector;}

	if(ncy==1){nsy[0]=1;nsy[1]=nsy[0];nsy[2]=nsy[0];}
		else{
				nsy[0] = NCY/nscy;						
				lcs    = NCY%nscy;									if(lcs>0){nsy[0]++;}
				nsy[1] = (NCY-shifty1)/nscy;				nsy[1]++;			
				lcs    = (NCY-shifty1)%nscy;				if(lcs>0){nsy[1]++;}
				nsy[2] = (NCY-shifty2)/nscy;				nsy[2]++;			
				lcs    = (NCY-shifty2)%nscy;				if(lcs>0){nsy[2]++;}}

	//////////////////////// nsz //////////////////////////
	if(sym==1){NCZ=ncz;}
		else{NCZ=(ncz/2)+sector;}

	if((ncz==1)||(shape==4)){nsz[0]=1;nsz[1]=nsz[0];nsz[2]=nsz[0];}
		else{
				nsz[0]=NCZ/nscz;
				lcs=NCZ%nscz;												if(lcs>0){nsz[0]++;}
				nsz[1] = (NCZ-shiftz1)/nscz;				nsz[1]++;			
				lcs  = (NCZ-shiftz1)%nscz;					if(lcs>0){nsz[1]++;}
				nsz[2] = (NCZ-shiftz2)/nscz;				nsz[2]++;			
				lcs  = (NCZ-shiftz2)%nscz;					if(lcs>0){nsz[2]++;}}

	xsall[0] = nsx[0]*nsy[0]*nsz[0];
	xsall[1] = nsx[1]*nsy[1]*nsz[1];	
	xsall[2] = nsx[2]*nsy[2]*nsz[2];

	if(sym==0){xsall[0]=1;xsall[1]=1;xsall[2]=1;}

	nsall  = xsall[0] + xsall[1] + xsall[2];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::nsectortoijk(int s, int n, int& i, int& j, int& k)
{
		// n = i + ncx*j + nx*ny*k;
		k=  n/(nsx[s]*nsy[s]);
		j= (n - k*nsx[s]*nsy[s])/nsx[s];
		i= (n - k*nsx[s]*nsy[s]) - (j*nsx[s]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Interpolation /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::interX(double r[3], int n)
{int i, j, k;
 double ax,xj,xjp1,xjm1,  ymax, ymax1, ymax2,  ymin, ymin1, ymin2,  zmax, zmax1, zmax2,  zmin, zmin1, zmin2;

	i=vsX[n].adr[0];
	j=vsX[n].adr[1];
	k=vsX[n].adr[2];

  ymax = vsY[ ijktonsurY(0,j+1,k) ].rc[1];	
	ymax1 = ymax - tola;
	ymax2 = ymax + tola;
	if(r[1]>ymax2) {return 0;}

	ymin = vsY[ ijktonsurY(0,j,k) ].rc[1];
	ymin1 = ymin - tola;
	ymin2 = ymin + tola;	
	if(r[1]<ymin1) {return 0;}

	zmax = vsZ[ ijktonsurZ(0,j,k+1) ].rc[2];
	zmax1 = zmax - tola;
	zmax2 = zmax + tola;
	if(r[2]>zmax2) {return 0;}

	zmin = vsZ[ ijktonsurZ(0,j,k) ].rc[2];
	zmin1 = zmin - tola;
	zmin2 = zmin + tola;
	if(r[2]<zmin1) {return 0;}

	xj=vsX[n].rc[0];

	if(i>0){xjm1=vsX[ijktonsurX(i-1,j,k)].rc[0];}	
			else{xjm1=-1.0;}

	if(i<ncx){xjp1=vsX[ijktonsurX(i+1,j,k)].rc[0];}	
			else{xjp1=x+1.0;}

	if((r[0]<=xjp1+tola) && (r[0]>=xjm1-tola))
		{if(r[0]>=xj) {ax=(xjp1-r[0])/(xjp1-xj);}
				else ax=(r[0]-xjm1)/(xj-xjm1);}
		else{return 0;}		

	if( (((r[1]>ymin1)&&(r[1]<ymin2))&&((r[2]>zmin1)&&(r[2]<zmin2))) ||
		  (((r[1]<ymax2)&&(r[1]>ymax1))&&((r[2]>zmin1)&&(r[2]<zmin2))) || 
			(((r[1]<ymax2)&&(r[1]>ymax1))&&((r[2]<zmax2)&&(r[2]>zmax1))) || 
			(((r[1]>ymin1)&&(r[1]<ymin2))&&((r[2]<zmax2)&&(r[2]>zmax1))) )
		{return ax*0.25;}

	if( (((r[1]>ymin1)&&(r[1]<ymin2)) || ((r[1]<ymax2)&&(r[1]>ymax1)) || ((r[2]>zmin1)&&(r[2]<zmin2)) || ((r[2]<zmax2)&&(r[2]>zmax1))) &&
			  (r[1]>tola) && (r[1]<(y-tola)) && (r[2]>tola) && (r[2]<(z-tola))  )	
			{return ax*0.5;}	

	return ax;	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::interY(double r[3], int n)
{int i, j, k;
 double ay,yj,yjp1,yjm1,  xmax, xmax1, xmax2,  xmin, xmin1, xmin2,  zmax, zmax1, zmax2,  zmin, zmin1, zmin2;

	i=vsY[n].adr[0];
	j=vsY[n].adr[1];
	k=vsY[n].adr[2];

	xmax = vsX[ ijktonsurX(i+1,0,k) ].rc[0];
	xmax1 = xmax - tola;
	xmax2 = xmax + tola;
	if(r[0]>xmax2) {return 0;}
	
	xmin = vsX[ ijktonsurX(i,0,k) ].rc[0];
	xmin1 = xmin - tola;
	xmin2 = xmin + tola;	
	if(r[0]<xmin1) {return 0;}

	zmax = vsZ[ ijktonsurZ(i,0,k+1) ].rc[2];
	zmax1 = zmax - tola;
	zmax2 = zmax + tola;
	if(r[2]>zmax2) {return 0;}

	zmin = vsZ[ ijktonsurZ(i,0,k) ].rc[2];
	zmin1 = zmin - tola;
	zmin2 = zmin + tola;
	if(r[2]<zmin1) {return 0;}

	yj=vsY[n].rc[1];

	if(j>0){yjm1=vsY[ijktonsurY(i,j-1,k)].rc[1];}	
			else{yjm1=-1.0;}

	if(j<ncy){yjp1=vsY[ijktonsurY(i,j+1,k)].rc[1];}	
			else{yjp1=y+1.0;}

	if ((r[1]<=yjp1+tola)&&(r[1]>=yjm1-tola))
		{if(r[1]>=yj) ay=(yjp1-r[1])/(yjp1-yj);
					else ay = (r[1]-yjm1)/(yj-yjm1);}
			else{ay=0;}		

	if	( (((r[0]>xmin1)&&(r[0]<xmin2))&&((r[2]>zmin1)&&(r[2]<zmin2))) ||
			 	(((r[0]<xmax2)&&(r[0]>xmax1))&&((r[2]>zmin1)&&(r[2]<zmin2))) || 
				(((r[0]<xmax2)&&(r[0]>xmax1))&&((r[2]<zmax2)&&(r[2]>zmax1))) || 
				(((r[0]>xmin1)&&(r[0]<xmin2))&&((r[2]<zmax2)&&(r[2]>zmax1))))
		{return ay*0.25;}
		
	if( ( ((r[0]>xmin1)&&(r[0]<xmin2)) || ((r[0]<xmax2)&&(r[0]>xmax1)) || ((r[2]>zmin1)&&(r[2]<zmin2)) || ((r[2]<zmax2)&&(r[2]>zmax1)) ) &&
			 (r[0]>=tola) && (r[0]<=(x-tola)) && (r[2]>=tola) && (r[2]<=(z-tola)) ) 
		{return ay*0.5;}

	return ay;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::interZ(double r[3], int n)
{int i, j, k;
 double az,zj,zjp1,zjm1,  ymax, ymax1, ymax2,  ymin, ymin1, ymin2,  xmax, xmax1, xmax2,  xmin, xmin1, xmin2;
	
	i=vsZ[n].adr[0];
	j=vsZ[n].adr[1];
	k=vsZ[n].adr[2];

	ymax = vsY[ ijktonsurY(i,j+1,0) ].rc[1];
	ymax1 = ymax - tola;
	ymax2 = ymax + tola;
	if(r[1]>ymax2) {return 0;}

	ymin = vsY[ ijktonsurY(i,j,0) ].rc[1];
	ymin1 = ymin - tola;
	ymin2 = ymin + tola;	
	if(r[1]<ymin1) {return 0;}

	xmax = vsX[ ijktonsurX(i+1,j,0) ].rc[0];
	xmax1 = xmax - tola;
	xmax2 = xmax + tola;
	if(r[0]>xmax2) {return 0;}

	xmin = vsX[ ijktonsurX(i,j,0) ].rc[0];
	xmin1 = xmin - tola;
	xmin2 = xmin + tola;	
	if(r[0]<xmin1) {return 0;}

	zj=vsZ[n].rc[2];

	if(k>0){zjm1=vsZ[ijktonsurZ(i,j,k-1)].rc[2];}	
			else{zjm1=-1.0;}

	if(k<ncz){zjp1=vsZ[ijktonsurZ(i,j,k+1)].rc[2];}	
			else{zjp1=z+1.0;}

	if( (r[2]<=zjp1+tola) && (r[2]>=zjm1-tola) )
		{if(r[2]>=zj) az=(zjp1-r[2])/(zjp1-zj);		
			else az=(r[2]-zjm1)/(zj-zjm1);}
			else{az=0;}		

	if( (((r[1]>ymin1)&&(r[1]<ymin2))&&((r[0]>xmin1)&&(r[0]<xmin2))) ||
			(((r[1]<ymax2)&&(r[1]>ymax1))&&((r[0]>xmin1)&&(r[0]<xmin2))) ||
			(((r[1]<ymax2)&&(r[1]>ymax1))&&((r[0]<xmax2)&&(r[0]>xmax1))) ||
			(((r[1]>ymin1)&&(r[1]<ymin2))&&((r[0]<xmax2)&&(r[0]>xmax1))))
		{	return az*0.25;}
	if( (((r[1]>ymin1)&&(r[1]<ymin2)) || ((r[1]<ymax2)&&(r[1]>ymax1)) || ((r[0]>xmin1)&&(r[0]<xmin2)) || ((r[0]<xmax2)&&(r[0]>xmax1))) && 
			(r[0]>tola) && (r[0]<(x-tola)) && (r[1]>tola) && (r[1]<(y-tola)) )
			{return az*0.5;}
	
	return az;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::interJ(double r[3], double J[3])
{int n=0;
 double jx=0, jy=0, jz=0, w[3]={0};

	if (ncx!=1){
			for(n=0;n<nsurx;n++){
									vec_subs2(vsX[n].rc,r,w);
									vec_fabs(w);

									if ((w[0]<=vsX[n].size[0] + tola)&&(w[1]<=vsX[n].size[1]/2.0 + tola)&&(w[2]<=vsX[n].size[2] + tola))
										 {jx+= interX(r,n) * (vsX[n].Jo + vsX[n].dJ);}}}
					else{jx=0;}

	if (ncy!=1){
			for(n=0;n<nsury;n++){
									vec_subs2(vsY[n].rc,r,w);
									vec_fabs(w);

									if ((w[0]<=vsY[n].size[0]/2.0 + tola)&&(w[1]<=vsY[n].size[1] + tola)&&(w[2]<=vsY[n].size[2] + tola))
										 {jy+= interY(r,n) * (vsY[n].Jo + vsY[n].dJ);}}}
					else{jy=0;}

	if (ncz!=1){
			for(n=0;n<nsurz;n++){
									vec_subs2(vsZ[n].rc,r,w);	
									vec_fabs(w);

									if ((w[0]<=vsZ[n].size[0]/2.0 + tola)&&(w[1]<=vsZ[n].size[1]/2.0 + tola)&&(w[2]<=vsZ[n].size[2] + tola))
										 {jz+= interZ(r,n) * (vsZ[n].Jo + vsZ[n].dJ);}}}
					else{jz=0;}	
		
	J[0]=jx; 
	J[1]=jy;
	J[2]=jz;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::interJ_cell(int n, double J[3])
{int i,j,k,n1,n2;

	i=vc[n].adr[0];
  j=vc[n].adr[1];
  k=vc[n].adr[2];

	n1=ijktonsurX(i,j,k);
	n2=ijktonsurX(i+1,j,k);
	J[0] = (1.0/2.0)*(vsX[n1].Jo + vsX[n1].dJ	+	vsX[n2].Jo + vsX[n2].dJ);

	J[1] = (1.0/2.0)*(vsY[ijktonsurY(i,j,k)].Jo + vsY[ijktonsurY(i,j,k)].dJ + vsY[ijktonsurY(i,j+1,k)].Jo + vsY[ijktonsurY(i,j+1,k)].dJ);
	if(Ismax!=0) {J[1]+= vc[n].Js[1];}

	J[2] = (1.0/2.0)*(vsZ[ijktonsurZ(i,j,k)].Jo + vsZ[ijktonsurZ(i,j,k)].dJ	+	vsZ[ijktonsurZ(i,j,k+1)].Jo + vsZ[ijktonsurZ(i,j,k+1)].dJ);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::interJ0_cell(int n, double J[3])
{int i,j,k;

	i=vc[n].adr[0];
  j=vc[n].adr[1];
  k=vc[n].adr[2];

	J[0] = (1.0/2.0)*(vsX[ijktonsurX(i,j,k)].Jo +	vsX[ijktonsurX(i+1,j,k)].Jo);
	J[1] = (1.0/2.0)*(vsY[ijktonsurY(i,j,k)].Jo + vsY[ijktonsurY(i,j+1,k)].Jo);
	J[2] = (1.0/2.0)*(vsZ[ijktonsurZ(i,j,k)].Jo +	vsZ[ijktonsurZ(i,j,k+1)].Jo);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::interJs_cell(int n, int s)
{int i,j,k,n1,n2;

	i=vc[n].adr[0];
  j=vc[n].adr[1];
  k=vc[n].adr[2];

	if(s==0) {return 0;}

	if(s==1){
			n1=ijktonsurY(i,j,k);
			n2=ijktonsurY(i,j+1,k);

			return (1.0/2.0)*(vsY[n1].Is*vsY[n1].us + vsY[n2].Is*vsY[n2].us);}

	if(s==2) {return 0;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::intercell_A(int n, double A[3])
{int i,j,k;

	i=vc[n].adr[0];
  j=vc[n].adr[1];
  k=vc[n].adr[2];

	A[0] = (1.0/2.0)*(vsX[ijktonsurX(i,j,k)].Aa0 + vsX[ijktonsurX(i,j,k)].dAa + vsX[ijktonsurX(i+1,j,k)].Aa0 + vsX[ijktonsurX(i+1,j,k)].dAa);
	A[1] = (1.0/2.0)*(vsY[ijktonsurY(i,j,k)].Aa0 + vsY[ijktonsurY(i,j,k)].dAa + vsY[ijktonsurY(i,j+1,k)].Aa0 + vsY[ijktonsurY(i,j+1,k)].dAa);
	A[2] = (1.0/2.0)*(vsZ[ijktonsurZ(i,j,k)].Aa0 + vsZ[ijktonsurZ(i,j,k)].dAa + vsZ[ijktonsurZ(i,j,k+1)].Aa0 + vsZ[ijktonsurZ(i,j,k+1)].dAa);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::intercell_T(int n, double T[3])
{int i,j,k;
 double dT[3]={0},T0[3]={0};	

	i=vc[n].adr[0];
  j=vc[n].adr[1];
  k=vc[n].adr[2];

	dT[0] = (1.0/4.0)*(veX[ijktonedgeX(i,j,k)].dT + veX[ijktonedgeX(i,j+1,k)].dT + veX[ijktonedgeX(i,j,k+1)].dT + veX[ijktonedgeX(i,j+1,k+1)].dT);
	dT[1] = (1.0/4.0)*(veY[ijktonedgeY(i,j,k)].dT + veY[ijktonedgeY(i+1,j,k)].dT + veY[ijktonedgeY(i,j,k+1)].dT + veY[ijktonedgeY(i+1,j,k+1)].dT);
	dT[2] = (1.0/4.0)*(veZ[ijktonedgeZ(i,j,k)].dT + veZ[ijktonedgeZ(i+1,j,k)].dT + veZ[ijktonedgeZ(i,j+1,k)].dT + veZ[ijktonedgeZ(i+1,j+1,k)].dT);

	T0[0] = (1.0/4.0)*(veX[ijktonedgeX(i,j,k)].T0 + veX[ijktonedgeX(i,j+1,k)].T0 + veX[ijktonedgeX(i,j,k+1)].T0 + veX[ijktonedgeX(i,j+1,k+1)].T0);
	T0[1] = (1.0/4.0)*(veY[ijktonedgeY(i,j,k)].T0 + veY[ijktonedgeY(i+1,j,k)].T0 + veY[ijktonedgeY(i,j,k+1)].T0 + veY[ijktonedgeY(i+1,j,k+1)].T0);
	T0[2] = (1.0/4.0)*(veZ[ijktonedgeZ(i,j,k)].T0 + veZ[ijktonedgeZ(i+1,j,k)].T0 + veZ[ijktonedgeZ(i,j+1,k)].T0 + veZ[ijktonedgeZ(i+1,j+1,k)].T0);

	vec_add2(T0, dT, T);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::interpolation_cell()
{int n=0;

	#pragma omp parallel for private(n) num_threads(num_threads)
	for (n=0;n<nc;n++){
			interJ_cell(n, vc[n].J); 
			magnetic_moment(n,vc[n].m);
			intercell_A(n, vc[n].A);
			intercell_T(n, vc[n].T);		
			local_loss(n);}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::polar_cart_X(double rc[3], double Xp[3], double X[3])
{
	X[0] = Xp[0]*cos(rc[1]) - Xp[1]*sin(rc[1]);
	X[1] = Xp[0]*sin(rc[1]) + Xp[1]*cos(rc[1]);
	X[2] = Xp[2];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::bilinear_inter_JcB(double B, double theta)
{int n, n1, n2, s1, s2;
 double max=10000,v,ht1,ht2,hb1,hb2,tol;

	tol=1e-15;

	if(measured_points==1){ht1=1;s1=0;ht2=0;s2=0;}else{
				for(n=0;n<measured_points-1;n++){
								n1=n; n2=n+1;

								v = fabs(data[0][n1][2] - theta);   if(v<tol) {ht1=1; ht2=0; s1=n1; s2=0; }
								v = fabs(data[0][n2][2] - theta);   if(v<tol) {ht1=0; ht2=1; s1=0;  s2=n2;}

								if((data[0][n1][2]<theta)&&(theta<data[0][n2][2])){
															c   = fabs(data[0][n1][2] - data[0][n2][2]);

															ht2 = fabs(data[0][n1][2] - theta)/c;
															ht1 = fabs(data[0][n2][2] - theta)/c;
															s1=n1; s2=n2;}}}

	for(n=0;n<measured_fields-1;n++){
						n1=n; n2=n+1;

						v = fabs(data[n1][0][1] - B);   if(v<tol) {return data[n1][s1][3]*ht1 + data[n1][s2][3]*ht2;}
						v = fabs(data[n2][0][1] - B);   if(v<tol) {return data[n2][s1][3]*ht1 + data[n2][s2][3]*ht2;}

						if((data[n1][s1][1]<B)&&(B<data[n2][s2][1])){
													c   = fabs(data[n1][0][1] - data[n2][0][1]);

													hb2 = fabs(data[n1][0][1] - B)/c;
													hb1 = fabs(data[n2][0][1] - B)/c;

					return (data[n1][s1][3]*ht1 + data[n1][s2][3]*ht2) * hb1  +  (data[n2][s1][3]*ht1 + data[n2][s2][3]*ht2) * hb2;}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::bilinear_inter_nB(double B, double theta)
{int n, n1, n2, s1, s2;
 double max=10000,v,ht1,ht2,hb1,hb2,tol;

	tol=1e-15;

	if(measured_points==1){ht1=1;s1=0;ht2=0;s2=0;}else{
				for(n=0;n<measured_points-1;n++){
								n1=n; n2=n+1;

								v = fabs(data[0][n1][2] - theta);   if(v<tol) {ht1=1; ht2=0; s1=n1; s2=0; }
								v = fabs(data[0][n2][2] - theta);   if(v<tol) {ht1=0; ht2=1; s1=0;  s2=n2;}

								if((data[0][n1][2]<theta)&&(theta<data[0][n2][2])){
															c   = fabs(data[0][n1][2] - data[0][n2][2]);

															ht2 = fabs(data[0][n1][2] - theta)/c;
															ht1 = fabs(data[0][n2][2] - theta)/c;
															s1=n1; s2=n2;}}}

	for(n=0;n<measured_fields-1;n++){
						n1=n; n2=n+1;

						v = fabs(data[n1][0][1] - B);   if(v<tol) {return data[n1][s1][4]*ht1 + data[n1][s2][4]*ht2;}
						v = fabs(data[n2][0][1] - B);   if(v<tol) {return data[n2][s1][4]*ht1 + data[n2][s2][4]*ht2;}

						if((data[n1][s1][1]<B)&&(B<data[n2][s2][1])){
													c   = fabs(data[n1][0][1] - data[n2][0][1]);

													hb2 = fabs(data[n1][0][1] - B)/c;
													hb1 = fabs(data[n2][0][1] - B)/c;

					return (data[n1][s1][4]*ht1 + data[n1][s2][4]*ht2) * hb1  +  (data[n2][s1][4]*ht1 + data[n2][s2][4]*ht2) * hb2;}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Save data to memory and file //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_nodes() 
{int m,i,j,k; 
 double RC[6]={0};

	for(m=0;m<nn;m++){												
					nnodetoijk(m, i, j, k);
		      vnode[m].adr[0]=i;
		      vnode[m].adr[1]=j;
		      vnode[m].adr[2]=k;       

					node_rc(i, j, k, RC);  
					vec_cpy_n(6,vnode[m].rc, RC);}

	fstream f;
	f.open("data/nodes.txt",ios::out);
	f << "#1	2	 3  4  5	6	 7	8	 9  10" << endl; 
	f << "#n  i  j  k	 x	y  z  r	 fi  z" << endl; 
	for(m=0;m<nn;m++){												
				f<< m <<"\t"<< vnode[m].adr[0] <<"\t"<< vnode[m].adr[1] <<"\t"<< vnode[m].adr[2] <<"\t"<< vnode[m].rc[0] <<"\t"<< vnode[m].rc[1] <<"\t"<< vnode[m].rc[2] <<
								"\t"<< vnode[m].rc[3] <<"\t"<< vnode[m].rc[4] <<"\t"<< vnode[m].rc[5] <<endl<<flush;}
	f.close();
}		

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_edges() 
{int m,i,j,k,x1,x2,y1,y2,z1,z2,m1,m2,m3,m4; 
 double RC[6]={0}, l=0;

	for(m=0;m<nedgex;m++){												
					nedgeXtoijk(m, i, j, k);
		      veX[m].adr[0] = i;
		      veX[m].adr[1] = j;
		      veX[m].adr[2] = k;    

					if((k!=0)||(j!=0)||(k!=ncz)||(j!=ncy)){
							edgeXtosurYZ(m,y1,y2,z1,z2);
							veX[m].adrs[0] = y1;
							veX[m].adrs[1] = y2;
							veX[m].adrs[2] = z1;
							veX[m].adrs[3] = z2;

							veX[m].adrc[0] = icell(i, j, k-1);
							veX[m].adrc[1] = icell(i,j,k);
							veX[m].adrc[2] = icell(i, j-1, k);
							veX[m].adrc[3] = icell(i,j-1,k-1);

							veX[m].adres[0] = ijktonedgeX(ncx-1-i,j,k);
							veX[m].adres[1] = ijktonedgeX(i,ncy-j,k);
							veX[m].adres[2] = ijktonedgeX(ncx-1-i,ncy-j,k);
							veX[m].adres[3] = ijktonedgeX(i,j,ncz-k);
							veX[m].adres[4] = ijktonedgeX(ncx-1-i,j,ncz-k);	
							veX[m].adres[5] = ijktonedgeX(i,ncy-j,ncz-k);	
							veX[m].adres[6] = ijktonedgeX(ncx-1-i,ncy-j,ncz-k);}

					edgeX_rc(i, j, k, RC);													
				 	vec_cpy_n(6,veX[m].rc, RC);
					veX[m].l = edgeX_lenght(i,j,k); 							
					veX[m].T0 = 0;																
					veX[m].dT = 0;		
					veX[m].dTx= 0;}

	for(m=0;m<nedgey;m++){												
					nedgeYtoijk(m, i, j, k);
		      veY[m].adr[0]=i;
		      veY[m].adr[1]=j;
		      veY[m].adr[2]=k;  

					if((i!=0)||(k!=0)||(i!=ncx)||(k!=ncz)){
							edgeYtosurXZ(m,x1,x2,z1,z2);
							veY[m].adrs[0] = x1;
							veY[m].adrs[1] = x2;
							veY[m].adrs[2] = z1;
							veY[m].adrs[3] = z2;

							veY[m].adrc[0] = icell(i, j, k-1);
							veY[m].adrc[1] = icell(i,j,k);
							veY[m].adrc[2] = icell(i-1,j, k);
							veY[m].adrc[3] = icell(i-1,j,k-1);

							veY[m].adres[0] = ijktonedgeY(ncx-i,j,k);
							veY[m].adres[1] = ijktonedgeY(i,ncy-1-j,k);
							veY[m].adres[2] = ijktonedgeY(ncx-i,ncy-1-j,k);
							veY[m].adres[3] = ijktonedgeY(i,j,ncz-k);
							veY[m].adres[4] = ijktonedgeY(ncx-i,j,ncz-k);
							veY[m].adres[5] = ijktonedgeY(i,ncy-1-j,ncz-k);
							veY[m].adres[6] = ijktonedgeY(ncx-i,ncy-1-j,ncz-k);}

					edgeY_rc(i, j, k, RC);													
				 	vec_cpy_n(6,veY[m].rc, RC);
					veY[m].l = edgeY_lenght(i,j,k); 							
					veY[m].T0 = 0;	
					veY[m].dT = 0;		
					veY[m].dTx= 0;}

	for(m=0;m<nedgez;m++){												
					nedgeZtoijk(m, i, j, k);
		      veZ[m].adr[0]=i;
		      veZ[m].adr[1]=j;
		      veZ[m].adr[2]=k; 

					if((i!=0)||(j!=0)||(i!=ncx)||(j!=ncy)){
							edgeZtosurXY(m,x1,x2,y1,y2);
							veZ[m].adrs[0] = x1;
							veZ[m].adrs[1] = x2;
							veZ[m].adrs[2] = y1;
							veZ[m].adrs[3] = y2;

							veZ[m].adrc[0] = icell(i, j-1, k);
							veZ[m].adrc[1] = icell(i,j,k);
							veZ[m].adrc[2] = icell(i-1, j,k);
							veZ[m].adrc[3] = icell(i-1,j-1,k);

							veZ[m].adres[0] = ijktonedgeZ(ncx-i,j,k);
							veZ[m].adres[1] = ijktonedgeZ(i,ncy-j,k);
							veZ[m].adres[2] = ijktonedgeZ(ncx-i,ncy-j,k);
							veZ[m].adres[3] = ijktonedgeZ(i,j,ncz-1-k);
							veZ[m].adres[4] = ijktonedgeZ(ncx-i,j,ncz-1-k);
							veZ[m].adres[5] = ijktonedgeZ(i,ncy-j,ncz-1-k);
							veZ[m].adres[6] = ijktonedgeZ(ncx-i,ncy-j,ncz-1-k);}

					edgeZ_rc(i, j, k, RC);											
				 	vec_cpy_n(6,veZ[m].rc, RC);
					veZ[m].l = edgeZ_lenght(i,j,k);  						
					veZ[m].T0 = 0;					
					veZ[m].dT = 0;			
					veZ[m].dTx= 0;}

	if(shape==5){boundary_addres_edge();}							//calculates adreesses of cells and surfaces around the edges at the boundary

	fstream f;
	f.open("data/edgesX.txt",ios::out|ios::app);
	f << "#1	2	 3	4	 5	6	 7	8  9  10" << endl; 	
	f << "#n  i  j  k  x  y  z  r  fi z" << endl; 

	for(m=0; m<nedgex; m++){		
			f<< m <<"\t"<< veX[m].adr[0] <<"\t"<< veX[m].adr[1] <<"\t"<< veX[m].adr[2] <<"\t"<< veX[m].rc[0] <<"\t"<< veX[m].rc[1] <<"\t"<< veX[m].rc[2]
					  <<"\t"<< veX[m].rc[3] <<"\t"<< veX[m].rc[4] <<"\t"<< veX[m].rc[5] << endl << flush;}
	f.close();


	f.open("data/edgesY.txt",ios::out|ios::app);
	f << "#1	2	 3	4	 5	6	 7	8  9  10" << endl; 	
	f << "#n  i  j  k  x  y  z  r  fi z" << endl; 

	for(m=0; m<nedgey; m++){												
			f<< m <<"\t"<< veY[m].adr[0] <<"\t"<< veY[m].adr[1] <<"\t"<< veY[m].adr[2] <<"\t"<< veY[m].rc[0] <<"\t"<< veY[m].rc[1] <<"\t"<< veY[m].rc[2]
						<<"\t"<< veY[m].rc[3] <<"\t"<< veY[m].rc[4] <<"\t"<< veY[m].rc[5] << endl << flush;}
	f.close();


	f.open("data/edgesZ.txt",ios::out|ios::app);
	f << "#1	2	 3	4	 5	6	 7	8  9  10" << endl; 	
	f << "#n  i  j  k  x  y  z  r  fi z" << endl; 

	for(m=0; m<nedgez; m++){												
			f<< m <<"\t"<< veZ[m].adr[0] <<"\t"<< veZ[m].adr[1] <<"\t"<< veZ[m].adr[2] <<"\t"<< veZ[m].rc[0] <<"\t"<< veZ[m].rc[1] <<"\t"<< veZ[m].rc[2]
						<<"\t"<< veZ[m].rc[3] <<"\t"<< veZ[m].rc[4] <<"\t"<< veZ[m].rc[5] << endl << flush;}
	f.close();
}	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_surfaces() 
{int n, i, j, k; 
 double RC[6]={0}, SS[4]={0};

	for(n=0; n<nsurx; n++){
					nsurXtoijk (n, i, j, k);
		      vsX[n].adr[0] = i;
		      vsX[n].adr[1] = j;
		      vsX[n].adr[2] = k;    

					vsX[n].adre[0] = ijktonedgeY(i,j,k);
					vsX[n].adre[1] = ijktonedgeY(i,j,k+1);
					vsX[n].adre[2] = ijktonedgeZ(i,j+1,k);
					vsX[n].adre[3] = ijktonedgeZ(i,j,k);

					vsX[n].adrss[0] = ijktonsurX(ncx-i,j,k);
					vsX[n].adrss[1] = ijktonsurX(i,ncy-1-j,k);
					vsX[n].adrss[2] = ijktonsurX(ncx-i,ncy-1-j,k);
					vsX[n].adrss[3] = ijktonsurX(i,j,ncz-1-k);
					vsX[n].adrss[4] = ijktonsurX(ncx-i,j,ncz-1-k);	
					vsX[n].adrss[5] = ijktonsurX(i,ncy-1-j,ncz-1-k);
					vsX[n].adrss[6] = ijktonsurX(ncx-i,ncy-1-j,ncz-1-k);
				

					vsX[n].Vi = vol_infX(i, j, k);					//rework!!!!	

					surfacex_rc (i, j, k, RC);								
					vec_cpy_n(6,vsX[n].rc, RC);

					surfaceX_size ( i, j, k, SS);
					vec_cpy_n(4,vsX[n].size, SS);			
				
					vsX[n].S = surX_S(n,i,j,k);
	
					vsX[n].Jo      = 0;
					vsX[n].dJ      = 0;	
					vsX[n].dJx     = 0;		
					vsX[n].dAa     = 0;	
					vsX[n].Aa0     = 0;
					vsX[n].av_dA   = 0;
					vsX[n].dAx     = 0;
					vsX[n].us   	 = 0;	
					vsX[n].Is   	 = 0;	
					vsX[n].dI   	 = 0;	
					vsX[n].dAs   	 = 0;}

	for(n=0; n<nsury; n++){
					nsurYtoijk (n, i, j, k);
		      vsY[n].adr[0] = i;
		      vsY[n].adr[1] = j;
		      vsY[n].adr[2] = k; 

					vsY[n].adre[0] = ijktonedgeX(i,j,k);
		 		  vsY[n].adre[1] = ijktonedgeX(i,j,k+1);
			 	  vsY[n].adre[2] = ijktonedgeZ(i+1,j,k);
				  vsY[n].adre[3] = ijktonedgeZ(i,j,k);

					vsY[n].adrss[0] = ijktonsurY(ncx-1-i,j,k);
		 		  vsY[n].adrss[1] = ijktonsurY(i,ncy-j,k);
			 	  vsY[n].adrss[2] = ijktonsurY(ncx-1-i,ncy-j,k);	
				  vsY[n].adrss[3] = ijktonsurY(i,j,ncz-1-k);
					vsY[n].adrss[4] = ijktonsurY(ncx-1-i,j,ncz-1-k);
		 		  vsY[n].adrss[5] = ijktonsurY(i,ncy-j,ncz-1-k);
			 	  vsY[n].adrss[6] = ijktonsurY(ncx-1-i,ncy-j,ncz-1-k);

					vsY[n].Vi = vol_infY(i, j, k);						//rework!!!!

					surfacey_rc (i, j, k, RC);								

					vec_cpy_n(6,vsY[n].rc, RC);

					surfaceY_size ( i, j, k, SS);
					vec_cpy_n(4,vsY[n].size, SS);
					vsY[n].S = surY_S(n,i,j,k);

					vsY[n].Jo      = 0;
					vsY[n].dJ      = 0;	
					vsY[n].dJx     = 0;						
					vsY[n].dAa     = 0;	
					vsY[n].Aa0     = 0;
					vsY[n].av_dA   = 0;
					vsY[n].dAx     = 0;
					vsY[n].us	     = 0;
					vsY[n].Is   	 = 0;	
					vsY[n].dI   	 = 0;			
					vsY[n].dAs   	 = 0;}
 
	for(n=0; n<nsurz; n++){
					nsurZtoijk (n, i, j, k); 			
		      vsZ[n].adr[0] = i;
		      vsZ[n].adr[1] = j;
		      vsZ[n].adr[2] = k; 

					vsZ[n].adre[0] = ijktonedgeX(i,j,k);
					vsZ[n].adre[1] = ijktonedgeX(i,j+1,k);
					vsZ[n].adre[2] = ijktonedgeY(i+1,j,k);
					vsZ[n].adre[3] = ijktonedgeY(i,j,k);

					vsZ[n].adrss[0] = ijktonsurZ(ncx-1-i,j,k);
					vsZ[n].adrss[1] = ijktonsurZ(i,ncy-1-j,k);
					vsZ[n].adrss[2] = ijktonsurZ(ncx-1-i,ncy-1-j,k);
					vsZ[n].adrss[3] = ijktonsurZ(i,j,ncz-k);
					vsZ[n].adrss[4] = ijktonsurZ(ncx-1-i,j,ncz-k);
					vsZ[n].adrss[5] = ijktonsurZ(i,ncy-1-j,ncz-k);
					vsZ[n].adrss[6] = ijktonsurZ(ncx-1-i,ncy-1-j,ncz-k);

					vsZ[n].Vi = vol_infZ(i, j, k);							//rework!!!!

					surfacez_rc (i, j, k, RC);								
					vec_cpy_n(6,vsZ[n].rc, RC);

					surfaceZ_size ( i, j, k, SS);
					vec_cpy_n(4,vsZ[n].size, SS);

					vsZ[n].S = surZ_S(n,i,j,k);

					vsZ[n].Jo      = 0;
					vsZ[n].dJ			 = 0;
					vsZ[n].dJx     = 0;	
					vsZ[n].dAa     = 0;	
					vsZ[n].Aa0     = 0;
					vsZ[n].av_dA   = 0;	
					vsZ[n].dAx     = 0;	
					vsZ[n].us	     = 0;	
					vsZ[n].Is   	 = 0;	
					vsZ[n].dI   	 = 0;
					vsZ[n].dAs   	 = 0;}

		surf_shape();																				//set Jol, Nl for linear material according geometry

		fstream f;
		f.open("data/surfaces_X.txt",ios::out);
		f << "#1  2  3  4  5  6  7	8  9   10" << endl;
		f << "#n  i  j  k  x  y  z  r  FI  z" << endl;

		for(n=0; n<nsurx; n++){
					f<< n <<"\t"<< vsX[n].adr[0] <<"\t"<< vsX[n].adr[1] <<"\t"<< vsX[n].adr[2] <<"\t"<< vsX[n].rc[0] <<"\t"<< vsX[n].rc[1] <<"\t"<< vsX[n].rc[2] 
								<<"\t"<< vsX[n].rc[3]  <<"\t"<< vsX[n].rc[4]  <<"\t"<< vsX[n].rc[5]  << endl << flush;}
		f.close();

		f.open("data/surfaces_Y.txt",ios::out);
		f << "#1  2  3  4  5  6  7	8  9   10" << endl;
		f << "#n  i  j  k  x  y  z  r  FI  z" << endl; 

		for(n=0; n<nsury; n++){
					f<< n <<"\t"<< vsY[n].adr[0] <<"\t"<< vsY[n].adr[1] <<"\t"<< vsY[n].adr[2] <<"\t"<< vsY[n].rc[0] <<"\t"<< vsY[n].rc[1] <<"\t"<< vsY[n].rc[2] 
								<<"\t"<< vsY[n].rc[3]  <<"\t"<< vsY[n].rc[4]  <<"\t"<< vsY[n].rc[5]  << endl << flush;}
		f.close();

		f.open("data/surfaces_Z.txt",ios::out);
		f << "#1  2  3  4  5  6  7	8  9   10" << endl;
		f << "#n  i  j  k  x  y  z  r  FI  z" << endl;  

		for(n=0; n<nsurz; n++){
					f<< n <<"\t"<< vsZ[n].adr[0] <<"\t"<< vsZ[n].adr[1] <<"\t"<< vsZ[n].adr[2] <<"\t"<< vsZ[n].rc[0] <<"\t"<< vsZ[n].rc[1] <<"\t"<< vsZ[n].rc[2]
								<<"\t"<< vsZ[n].rc[3]  <<"\t"<< vsZ[n].rc[4]  <<"\t"<< vsZ[n].rc[5]  << endl << flush;}
		f.close();
}	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_cells() 
{int i,j, k, m;
 double SC[3]={0}, RC[6]={0};

	for(m=0; m<nc; m++){												
		      ncelltoijk (m, i, j, k);
		      vc[m].adr[0]=i;
		      vc[m].adr[1]=j;
		      vc[m].adr[2]=k;    

					cell_size (i, j, k, SC);
					vec_cpy(vc[m].size, SC);	
			
					vc[m].V = cell_vol (SC, i ,j ,k);
					cell_rc (i, j, k, RC);										
				  vec_cpy_n(6,vc[m].rc, RC);								

					vec_cpy(vc[m].J, zero);			 
					vec_cpy(vc[m].uv,zero);       
					vc[m].loss								=0;
					vc[m].N                   =N;
					vc[m].Jc		              =Jo;
					vc[m].Jcx		              =0;
					vc[m].U										=0;
					vec_cpy(vc[m].B_J, zero);
					vec_cpy(vc[m].A, zero);
					vec_cpy(vc[m].B, zero);
					vec_cpy(vc[m].pB, zero);
					vec_cpy(vc[m].T, zero);
					vec_cpy(vc[m].m, zero);
					vec_cpy(vc[m].E, zero);}

	if(Btrape==1){
			for(m=0;m<nc_plane;m++){
						ncellplanetoijk(m, i, j, k);

						RC_plane[m][0] = cell_rc_plane(i,j,k,0);
						RC_plane[m][1] = cell_rc_plane(i,j,k,1);
						RC_plane[m][2] = cell_rc_plane(i,j,k,2);

						vec_cpy(BJ_plane[m],zero);
						vec_cpy(B_plane[m],zero);}}

	if ((nncx!=0)||(shape==4)||(shape==3)) {cell_lin();}		//set Jol, Nl for linear material

	fstream f;
	f.open("data/cells.txt",ios::out);
	f << "#1	2	 3	4	 5	6	 7  8  9  10 11 12 13"  << endl; 	
	f << "#n  i	 j  k  x	y	 z	r	 FI z  i  j  k "  << endl; 

	for(m=0; m<nc; m++){	
			f<< m <<"\t"<< vc[m].adr[0] <<"\t"<< vc[m].adr[1] <<"\t"<< vc[m].adr[2] <<"\t"<< vc[m].rc[0] <<"\t"<< vc[m].rc[1] <<"\t"<< vc[m].rc[2] <<"\t"<< 
										 vc[m].rc[3] <<"\t"<< vc[m].rc[4] <<"\t"<< vc[m].rc[5] <<"\t"<< vc[m].adrp[0] <<"\t"<< vc[m].adrp[1] <<"\t"<< vc[m].adrp[2] <<endl<<flush;}	
	f.close();
}		

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Vector potential //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::surX_dA(int i)
{double ex[]={1.0,0,0}, xi=0, angle=0, Aa[3]={0};

	xi = vsX[i].rc[2] - 0.5*z; 
	angle = sin(fi*pi/180.0)*sin(theta*pi/180.0);

	vec_scale(ex, Ba);
	vec_scale(ex, xi);
	vec_scale(ex, angle);
	vec_cpy(Aa, ex); 
	return Aa[0] - vsX[i].Aa0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::surY_dA(int i)
{double ey[]={0,1.0,0},yi=0, angle=0, Aa[3]={0};

	yi=vsY[i].rc[0] - 0.5*x;

	angle = cos(fi*pi/180.0);

	vec_scale(ey, Ba);
	vec_scale(ey, yi);
	vec_scale(ey, angle);
	vec_cpy(Aa, ey); 

	return Aa[1] - vsY[i].Aa0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::surZ_dA(int i)
{double ez[]={0,0,1.0}, zi=0, angle=0, Aa[3]={0};

	zi=vsZ[i].rc[1] - 0.5*y;
	angle = sin(fi*pi/180.0)*cos(theta*pi/180.0);

	vec_scale(ez, Ba);
	vec_scale(ez, zi);
	vec_scale(ez, angle);
	vec_cpy(Aa, ez);

	return Aa[2] - vsZ[i].Aa0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
void class_cube::cal_av_dA(int s) 
{int m,n,n1,n2,i,j,k,o,p,r,l;
 double a;

	if(sym==1){
				if(ncx!=1){
								#pragma omp parallel for private(n,i,j,a) num_threads(num_threads)
								for(n=0;n<nsurx;n++){
												a=0.0;					
												for(i=0;i<nsurx;i++)	{a+=vsX[i].dJ * vsX[i].Vi * read_cmAx(i,n);}
												vsX[n].av_dA = a;}}

				if(ncy!=1){
								#pragma omp parallel for private(n,i,j,a) num_threads(num_threads)
								for(n=0;n<nsury;n++){
												a=0.0;
												for(i=0;i<nsury;i++)	{a+=vsY[i].dJ * vsY[i].Vi * read_cmAy(i,n);}			
												vsY[n].av_dA = a;}}

				if(ncz!=1){
								#pragma omp parallel for private(n,i,a) num_threads(num_threads)
								for(n=0;n<nsurz;n++){
												a=0.0;
												for(i=0;i<nsurz;i++)	{a+=vsZ[i].dJ * vsZ[i].Vi * read_cmAz(i,n);}
												vsZ[n].av_dA = a;}}}

	if(sym==2){
						s=s+1;
						if(s==3){s=0;}

						nrange(s);
						l=xsall[s]-1;

						if(ncx!=1){
										#pragma omp parallel for private(i) num_threads(num_threads)
										for(i=0;i<nsurx;i++) {vsX[i].av_dA = 0.0;} 
										
										#pragma omp parallel for if(ncz!=1) private(m,n,i,j,k,a) num_threads(num_threads)
										for(k=0;k<=cak[l];k++){													
											#pragma omp parallel for if(ncz==1) private(m,n,i,j,a) num_threads(num_threads)
											for(j=0;j<=caj[l];j++){
												for(i=0;i<=cai[l]+1;i++){
																	n = ijktonsurX(i,j,k);
																	a=0.0;					
																	for(m=0;m<nsurx;m++)	{a+=vsX[m].dJ * vsX[m].Vi * read_cmAx(m,n);}
																	vsX[n].av_dA = a;}}}}

						if(ncy!=1){
										#pragma omp parallel for private(i) num_threads(num_threads)
										for(i=0;i<nsury;i++) {vsY[i].av_dA = 0.0;}

										#pragma omp parallel for if(ncz!=1) private(m,n,i,j,k,a) num_threads(num_threads)
										for(k=0;k<=cak[l];k++){
											#pragma omp parallel for if(ncz==1) private(m,n,i,j,a) num_threads(num_threads)	
											for(j=0;j<=caj[l]+1;j++){
												for(i=0;i<=cai[l];i++){
																	n = ijktonsurY(i,j,k);
																	a=0.0;
																	for(m=0;m<nsury;m++)	{a+=vsY[m].dJ * vsY[m].Vi * read_cmAy(m,n);}				
																	vsY[n].av_dA = a;}}}}

						if(ncz!=1){
										#pragma omp parallel for private(i) num_threads(num_threads)
										for(i=0;i<nsurz;i++) {vsZ[i].av_dA = 0.0;}

										#pragma omp parallel for private(m,n,i,j,k,a) num_threads(num_threads)
										for(k=0;k<=cak[l]+1;k++){		
											for(j=0;j<=caj[l];j++){
												for(i=0;i<=cai[l];i++){
																	n = ijktonsurZ(i,j,k);
																	a=0.0;
																	for(m=0;m<nsurz;m++)	{a+=vsZ[m].dJ * vsZ[m].Vi * read_cmAz(m,n);}			
																	vsZ[n].av_dA = a;}}}}}
}	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::fast_matrix_cmA()
{int n;
 double R,x1,x2,x3,w[3]={0};

	if((cx>=cy)&&(cx>=cz))	{R=cx;}
	if((cy>=cx)&&(cy>=cz))	{R=cy;}
	if((cz>=cx)&&(cz>=cy))	{R=cz;}

	cout << "The analytical formula radius " << R << " m" << endl;

	cout << "cmAx..." << endl;
	#pragma omp parallel for private(n,x1,x2,x3) num_threads(num_threads)	 
	for(n=0;n<nsurx;n++){
		if(elc==0){
								x1=fabs(vsX[n].rc[0]-vsX[0].rc[0]);
								x2=fabs(vsX[n].rc[1]-vsX[0].rc[1]);
								x3=fabs(vsX[n].rc[2]-vsX[0].rc[2]);

								if((x1<R)&&(x2<R)&&(x3<R)){cmAx[n] = (ep0*mi0) * spvv(vsX[n].rc, vsX[0].rc, 1, vsX[n].size, vsX[0].size);}
										else{cmAx[n] = aprox_a(vsX[n].rc, vsX[0].rc);}																//approximation
							}else{cmAx[n] = num_a_auto(n,0,0);}}																								//elongated cells

	cout << "cmAy..." << endl;
	#pragma omp parallel for private(n,x1,x2,x3) num_threads(num_threads)	 
	for(n=0;n<nsury;n++){
		if(elc==0){
								x1=fabs(vsY[n].rc[0]-vsY[0].rc[0]);
								x2=fabs(vsY[n].rc[1]-vsY[0].rc[1]);
								x3=fabs(vsY[n].rc[2]-vsY[0].rc[2]);

								if((x1<R)&&(x2<R)&&(x3<R)){cmAy[n] = (ep0*mi0) * spvv(vsY[n].rc, vsY[0].rc, 1, vsY[n].size, vsY[0].size);}
										else{cmAy[n] = aprox_a(vsY[n].rc, vsY[0].rc);}															//approximation
							}else{cmAy[n] = num_a_auto(n,0,1);}}																							//elongated cells					

//		if((shape==5)&&(Bamax!=0)){cout <<"boundary cmAx,y " <<endl; boundary_cmA();}

	if(ncz!=1){
		cout << "cmAz..." << endl;
		#pragma omp parallel for private(n,x1,x2,x3) num_threads(num_threads)	 
			for(n=0;n<nsurz;n++){
				if(elc==0){
										x1=fabs(vsZ[n].rc[0]-vsZ[0].rc[0]);
										x2=fabs(vsZ[n].rc[1]-vsZ[0].rc[1]);
										x3=fabs(vsZ[n].rc[2]-vsZ[0].rc[2]);
										if((x1<R)&&(x2<R)&&(x3<R)){cmAz[n] = (ep0*mi0) * spvv(vsZ[n].rc, vsZ[0].rc, 1, vsZ[n].size, vsZ[0].size);}
												 else{cmAz[n] = aprox_a(vsZ[n].rc, vsZ[0].rc);}														//approximation
									}else{cmAz[n] = num_a_auto(n,0,2);}}}																						//elongated cells	

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::full_matrix_mA()
{int n,s,i,j,k,j1;
 double R,x1,x2,x3;

	s=1;
	
	if((vc[0].size[0]>=vc[0].size[1])&&(vc[0].size[0]>=vc[0].size[2]))	{R=vc[0].size[0]+vc[0].size[0]/300;}
	if((vc[0].size[1]>=vc[0].size[0])&&(vc[0].size[1]>=vc[0].size[2]))	{R=vc[0].size[1]+vc[0].size[1]/300;}
	if((vc[0].size[2]>=vc[0].size[0])&&(vc[0].size[2]>=vc[0].size[1]))	{R=vc[0].size[2]+vc[0].size[2]/300;}

	cout << "The analytical formula radius " << R << " m" << endl;
	cout << "mAx..." << endl;

	#pragma omp parallel for private(n,i,j,j1,x1,x2,x3) num_threads(num_threads)	 
	for(n=0;n<nsxall;n++){
			i=n/nsurx;
			j=n%nsurx;

			if(i<=j){
						x1=fabs(vsX[i].rc[0]-vsX[j].rc[0]);
						x2=fabs(vsX[i].rc[1]-vsX[j].rc[1]);
						x3=fabs(vsX[i].rc[2]-vsX[j].rc[2]);	

						j1=j-i;

						if(elc==0){
										if((x1<R)&&(x2<R)&&(x3<R))	
											{mAx[i][j1] = (ep0*mi0) * spvv(vsX[j].rc, vsX[i].rc,1,vsX[j].size, vsX[i].size);}
											else{mAx[i][j1] = aprox_a(vsX[i].rc, vsX[j].rc);}	

											}else{if(i==j){mAx[i][j1] = (ep0*mi0) * spvv(vsX[j].rc, vsX[i].rc,1,vsX[j].size, vsX[i].size);}
																else{mAx[i][j1] = num_a_auto(i,j,0);}}}}

	cout << "mAy..." << endl;
	#pragma omp parallel for private(n,i,j,j1,x1,x2,x3) num_threads(num_threads)	
	for(n=0;n<nsyall;n++){
			i=n/nsury;
			j=n%nsury;

			if(i<=j){
						x1=fabs(vsY[i].rc[0]-vsY[j].rc[0]);
						x2=fabs(vsY[i].rc[1]-vsY[j].rc[1]);
						x3=fabs(vsY[i].rc[2]-vsY[j].rc[2]);

						j1=j-i;

						if(elc==0){

										if((x1<R)&&(x2<R)&&(x3<R))	
											{mAy[i][j1] = (ep0*mi0) * spvv(vsY[j].rc, vsY[i].rc,1,vsY[j].size, vsY[i].size);}
											else{mAy[i][j1] = aprox_a(vsY[i].rc, vsY[j].rc);}

								}else{mAy[i][j1] = num_a_auto(i,j,1);}}}

	if(ncz!=1){
	cout << "mAz..." << endl;
	#pragma omp parallel for private(n,i,j,j1,x1,x2,x3) num_threads(num_threads)	
	for(n=0;n<nszall;n++){
			i=n/nsurz;
			j=n%nsurz;

			if(i<=j){
						x1=fabs(vsZ[i].rc[0]-vsZ[j].rc[0]);
						x2=fabs(vsZ[i].rc[1]-vsZ[j].rc[1]);
						x3=fabs(vsZ[i].rc[2]-vsZ[j].rc[2]);

						j1=j-i;

						if(elc==0){
										if((x1<R)&&(x2<R)&&(x3<R))  
											{mAz[i][j1] = (ep0*mi0) * spvv(vsZ[j].rc, vsZ[i].rc,1,vsZ[j].size, vsZ[i].size);}
											else{mAz[i][j1] = aprox_a(vsZ[i].rc, vsZ[j].rc);}	

								}else{mAz[i][j1] = num_a_auto(i,j,2);}}}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::read_cmAx(int i, int j)
{int i1=0,j1=0,l[3]={0},m[3]={0};

// i=n/nsurx;
// j=n%nsurx;

	if(full_matrix==0){
			vec_cpy_i(l,vsX[i].adr);
			vec_cpy_i(m,vsX[j].adr);

			return cmAx[ijktonsurX(abs(l[0]-m[0]),abs(l[1]-m[1]),abs(l[2]-m[2]))];}

	if(full_matrix==1){if(i<=j){return mAx[i][j-i];}else{return mAx[j][i-j];}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::read_cmAy(int i, int j)
{int i1=0,j1=0,k1=0,l[3]={0},m[3]={0};

// i=n/nsury;
// j=n%nsury;

	if(full_matrix==0){
			vec_cpy_i(l,vsY[i].adr);
			vec_cpy_i(m,vsY[j].adr);

			return cmAy[ijktonsurY(abs(l[0]-m[0]),abs(l[1]-m[1]),abs(l[2]-m[2]))];}
	
	if(full_matrix==1){if(i<=j){return mAy[i][j-i];}else{return mAy[j][i-j];}}	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::read_cmAz(int i, int j)
{int i1=0,j1=0,k1=0,l[3]={0},m[3]={0};

// i=n/nsurz;
// j=n%nsurz;

	if(full_matrix==0){
			vec_cpy_i(l,vsZ[i].adr);
			vec_cpy_i(m,vsZ[j].adr);

			return cmAz[ijktonsurZ(abs(l[0]-m[0]),abs(l[1]-m[1]),abs(l[2]-m[2]))];}
	
	if(full_matrix==1){if(i<=j){return mAz[i][j-i];}else{return mAz[j][i-j];}}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::num_a_auto(int l, int m, int i)
{int s=4,j, n[3]={3},na[3]={0}, nb[3]={0};
 double x1=0, x2, dx, dy, dz, a[3]={0}, b[3]={0}, c[3]={0};
	
	if(i==0)	{
			vec_cpy(b,vsX[l].size);
			vec_cpy(c,vsX[m].size);}
	if(i==1)	{
			vec_cpy(b,vsY[l].size);
			vec_cpy(c,vsY[m].size);}
	if(i==2)	{
			vec_cpy(b,vsZ[l].size);
			vec_cpy(c,vsZ[m].size);}
	do
		{
			if(fabs(l-m)<3){s*=2.0;}else{s++;}

			for(j=0;j<2;j++){
					if(j==0){vec_cpy(a,b);}else{vec_cpy(a,c);}

					if((a[0]>=a[1])&&(a[0]>=a[2])){
							dx 	 = a[0]/s;			
							n[0] = s;
							n[1] = int(a[1]/dx + 0.5); if(n[1]==0){n[1]=1;}	
							if(ncz==1){n[2]=1;}
								else{n[2] = int(a[2]/dx + 0.5); if(n[2]==0){n[2]=1;}}}

					if((a[1]>=a[0])&&(a[1]>=a[2])){
							dy	 = a[1]/s;			
							n[1] = s;
							n[0] = int(a[0]/dy + 0.5); if(n[0]==0){n[0]=1;}	
							if(ncz==1){n[2]=1;}
								else{n[2] = int(a[2]/dy + 0.5); if(n[2]==0){n[2]=1;}}}

					if((a[2]>=a[0])&&(a[2]>=a[1])){
							dz   = a[2]/s;			
							n[2] = s;
							n[0] = int(a[0]/dz + 0.5); if(n[0]==0){n[0]=1;}	
							n[1] = int(a[1]/dz + 0.5); if(n[1]==0){n[1]=1;}}

					if(j==0){vec_cpy_i(na,n);}else{vec_cpy_i(nb,n);}}

			x2=x1;
			x1=num_a(l, m, na, nb, i);
		}
		while((fabs(x1-x2)>x2*tol_elc)&&(na[0]<40)&&(na[1]<40));	

	return x1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::num_a(int l, int m, int *na, int *nb, int s)			
{int i[3]={0},j[3]={0},A[3]={1,1,1},B, adl[3]={0}, adm[3]={0};
 double rpp[3]={0}, x[3]={0}, x1[3]={0}, x0[3]={0}, x01[3]={0}, d[3]={0}, d1[3]={0}, Ax=0, Vs, Vs1, V1, V2, hi, hj, a[3]={0}, b[3]={0}, rcl[3]={0}, rcm[3]={0};

	if(s==0){
	vec_cpy(a,vsX[l].size);															//size of surface a(x,y,z)
	vec_cpy(b,vsX[m].size);															//size of surface b(x,y,z)
	vec_cpy(rcl,vsX[l].rc);															//position vector of surface l
	vec_cpy(rcm,vsX[m].rc);															//position vector of surface m
	vec_cpy_i(adl,vsX[l].adr);													//address of surface l
	vec_cpy_i(adm,vsX[m].adr);}													//address of surface m

	if(s==1){
	vec_cpy(a,vsY[l].size);															//size of surface a(x,y,z)
	vec_cpy(b,vsY[m].size);															//size of surface b(x,y,z)
	vec_cpy(rcl,vsY[l].rc);															//position vector of surface l
	vec_cpy(rcm,vsY[m].rc);															//position vector of surface m
	vec_cpy_i(adl,vsY[l].adr);													//address of surface l
	vec_cpy_i(adm,vsY[m].adr);}													//address of surface m

	if(s==2){
	vec_cpy(a,vsZ[l].size);															//size of surface a(x,y,z)
	vec_cpy(b,vsZ[m].size);															//size of surface b(x,y,z)
	vec_cpy(rcl,vsZ[l].rc);															//position vector of surface l
	vec_cpy(rcm,vsZ[m].rc);															//position vector of surface m
	vec_cpy_i(adl,vsZ[l].adr);													//address of surface l
	vec_cpy_i(adm,vsZ[m].adr);}													//address of surface m

	vec_div(d,a,na);																		//size of sub-element, na,nb-number of sub-elements
	vec_div(d1,b,nb);

	oper2(x0, rcl, a, d, s);														//zero position of sub-element
	oper2(x01, rcm, b, d1, s);																

	Vs = vol(d);																				//volume of sub-element
	Vs1= vol(d1);
	V1 = vol(a);																				//volume of influence of surface
	V2 = vol(b);

	if(s==0){
		if((fabs(adl[s]-adm[s])==1)&&(adl[1]==adm[1])&&(adl[2]==adm[2])) {A[s]=3;if(adl[s]<adm[s]) {vec_cpy(x01,x0);}else{vec_cpy(x0,x01);}}else{A[s]=2;}
		if((fabs(adl[s]-adm[s])<2)&&(adl[1]==adm[1])&&(adl[2]==adm[2])) {B=1;} else {B=0;}}

	if(s==1){
		if((fabs(adl[s]-adm[s])==1)&&(adl[0]==adm[0])&&(adl[2]==adm[2])) {A[s]=3;if(adl[s]<adm[s]) {vec_cpy(x01,x0);}else{vec_cpy(x0,x01);}}else{A[s]=2;}
		if((fabs(adl[s]-adm[s])<2)&&(adl[0]==adm[0])&&(adl[2]==adm[2])) {B=1;} else {B=0;}}

	if(s==2){
		if((fabs(adl[s]-adm[s])==1)&&(adl[0]==adm[0])&&(adl[1]==adm[1])) {A[s]=3;if(adl[s]<adm[s]) {vec_cpy(x01,x0);}else{vec_cpy(x0,x01);}}else{A[s]=2;}
		if((fabs(adl[s]-adm[s])<2)&&(adl[0]==adm[0])&&(adl[1]==adm[1])) {B=1;} else {B=0;}}

	for(i[0]=0;i[0]<na[0]*A[0];i[0]++){
			for(i[1]=0;i[1]<na[1]*A[1];i[1]++){
					for(i[2]=0;i[2]<na[2]*A[2];i[2]++){
								for(j[0]=0;j[0]<nb[0]*A[0];j[0]++){
										for(j[1]=0;j[1]<nb[1]*A[1];j[1]++){
												for(j[2]=0;j[2]<nb[2]*A[2];j[2]++){

																oper1(x, x0, i, d);																					//x0 + i*d change of sub-element position
																oper1(x1, x01, j, d1);

																hj=(-fabs(x[s] - rcl[s]) + a[s])/a[s]; if(hj<0){hj=0;}			//calculates h function
																hi=(-fabs(x1[s] - rcm[s]) + b[s])/b[s]; if(hi<0){hi=0;}

																vec_subs2(x,x1,rpp);																				//calculates distance of sub-elements

																if ((i[0]==j[0])&&(i[1]==j[1])&&(i[2]==j[2])&&(B==1))	{Ax+=4.0*pi*ep0*Vs*Vs1*hi*hj*spvv(x, x1, 1, d, d1);}	
																		else{Ax+=hj*hi*Vs*Vs1/vec_mod(rpp);}																		
																	}}}}}}	
	return mi0*Ax/(4.0*pi*V1*V2);		
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::matrix_snA()
{int i,n,m,s,m1[3]={0},m2[3]={0};
 double R;

	s=200;																														//50<=s<200

	cout << "Matrix for current source ... " << endl;

	if(ncx!=1){cout << "snAx..." << endl;
	#pragma omp parallel for private(n,i) num_threads(num_threads)	 
	for(n=0;n<nsurx;n++){
			snAx[n]=0;
			for(i=0;i<nsurx;i++){snAx[n]+= vsX[i].Vi * vsX[i].us * read_cmAx(i,n);}}}

	if(ncy!=1){cout << "snAy..." << endl;
	#pragma omp parallel for private(n,i,m) num_threads(num_threads)	 
	for(n=0;n<nsury;n++){
					snAy[n]=0;

					for(i=0;i<nsury;i++){snAy[n]+= vsY[i].Vi * vsY[i].us * read_cmAy(i,n);}

					if(shape==5){
									for(i=0;i<nsury;i++){															
													for(m=1;m<=s;m++)
																{
																	  {snAy[n]+= vsY[i].Vi * vsY[i].us * mi0/(4*pi*sqrt
																												(pow(vsY[n].rc[0]-vsY[i].rc[0],2) + pow(vsY[n].rc[1]-(vsY[i].rc[1]+(m*(y+cy))),2) + pow(vsY[n].rc[2]-vsY[i].rc[2],2)));}																
																															
																	  {snAy[n]+= vsY[i].Vi * vsY[i].us * mi0/(4*pi*sqrt
																												(pow(vsY[n].rc[0]-vsY[i].rc[0],2) + pow(vsY[n].rc[1]-(vsY[i].rc[1]-(m*(y+cy))),2) + pow(vsY[n].rc[2]-vsY[i].rc[2],2)));}	
																		}}}
						}}

	if(ncz!=1){	cout << "snAz..." << endl;
	#pragma omp parallel for private(n,i) num_threads(num_threads)	 
	for(n=0;n<nsurz;n++){
			snAz[n]=0;
			for(i=0;i<nsurz;i++){snAz[n]+= vsZ[i].Vi * vsZ[i].us * read_cmAz(i,n);}}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::boundary_cmA()
{int i,j,n,m,s;

	s=50;

	#pragma omp parallel for private(i,j,m) num_threads(num_threads)	 
	for(i=0;i<nsurx;i++){
			for(j=0;j<nsurx;j++){
						for(m=1;m<=s;m++){
										cmAx[i]+= mi0/(4*pi*sqrt( pow(vsX[i].rc[0]-vsX[j].rc[0],2)  +  pow(vsX[i].rc[1]-(vsX[j].rc[1]+(m*(y+cy))),2)  +  pow(vsX[i].rc[2]-vsY[j].rc[2],2) ));

										cmAx[i]+= mi0/(4*pi*sqrt( pow(vsX[i].rc[0]-vsX[j].rc[0],2)  +  pow(vsX[i].rc[1]-(vsX[j].rc[1]-(m*(y+cy))),2)  +  pow(vsX[i].rc[2]-vsY[j].rc[2],2) ));}}}



	#pragma omp parallel for private(i,j,m) num_threads(num_threads)	 
	for(i=0;i<nsury;i++){
			for(j=0;j<nsury;j++){
						for(m=1;m<=s;m++){
										cmAy[i]+= mi0/(4*pi*sqrt( pow(vsY[i].rc[0]-vsY[j].rc[0],2)  +  pow(vsY[i].rc[1]-(vsY[j].rc[1]+(m*(y+cy))),2)  +  pow(vsY[i].rc[2]-vsY[j].rc[2],2) ));

										cmAy[i]+= mi0/(4*pi*sqrt( pow(vsY[i].rc[0]-vsY[j].rc[0],2)  +  pow(vsY[i].rc[1]-(vsY[j].rc[1]-(m*(y+cy))),2)  +  pow(vsY[i].rc[2]-vsY[j].rc[2],2) ));}}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Magnetic field  ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::magnetic_field()
{int i,j;
 double Jxyz=0, Jxzy=0,  	Jyzx=0, Jyxz=0, Jzyz=0, 	Jzxy=0, Jzyx=0;

		#pragma omp parallel for private(i) num_threads(num_threads)	  
		for (i=0;i<nsurx;i++){vsX[i].J = vsX[i].Jo + vsX[i].dJ;}

		#pragma omp parallel for private(i) num_threads(num_threads)	  
		for (i=0;i<nsury;i++){vsY[i].J = vsY[i].Jo + vsY[i].dJ;}

		#pragma omp parallel for private(i) num_threads(num_threads)	  
		for (i=0;i<nsurz;i++){vsZ[i].J = vsZ[i].Jo + vsZ[i].dJ;}

		#pragma omp parallel for private(i,j,Jxyz,Jxzy,Jyxz,Jyzx,Jzxy,Jzyx) num_threads(num_threads)	  
		for (j=0;j<nc;j++){ 

						Jxyz=0, Jxzy=0, Jyzx=0, Jyxz=0, Jzyx=0, Jzxy=0;
						if((ncx!=1)||(cube==1)){
								for (i=0;i<nsurx;i++){
										Jyxz+= vsX[i].J * read_cmHyxz(j,i); 
										Jzxy+= vsX[i].J * read_cmHzxy(j,i);}}

						if((ncy!=1)||(cube==1)){
								for (i=0;i<nsury;i++){
										Jxyz+= vsY[i].J * read_cmHxyz(j,i); 
										Jzyx+= vsY[i].J * read_cmHzyx(j,i);}}

						if((ncz!=1)||(cube==1)){
								for (i=0;i<nsurz;i++){
										Jxzy+= vsZ[i].J * read_cmHxzy(j,i); 
										Jyzx+= vsZ[i].J * read_cmHyzx(j,i);}}

						vc[j].B_J[0] = (Jxyz - Jxzy); 																								
						vc[j].B[0] = vc[j].B_J[0] + B_a[0];
						vc[j].B_J[1] = (Jyzx - Jyxz);
						vc[j].B[1] = vc[j].B_J[1] + B_a[1];
						vc[j].B_J[2] = (Jzxy - Jzyx); 
						vc[j].B[2] = vc[j].B_J[2] + B_a[2];}	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::external_magnetic_field()
{int i,j;
 double Jxyz=0, Jxzy=0,		Jyzx=0, Jyxz=0,		Jzxy=0, Jzyx=0;

		#pragma omp parallel for private(i) num_threads(num_threads)	  
		for (i=0;i<nsurx;i++){vsX[i].J = vsX[i].Jo + vsX[i].dJ;}

		#pragma omp parallel for private(i) num_threads(num_threads)	  
		for (i=0;i<nsury;i++){vsY[i].J = vsY[i].Jo + vsY[i].dJ;}

		#pragma omp parallel for private(i) num_threads(num_threads)	  
		for (i=0;i<nsurz;i++){vsZ[i].J = vsZ[i].Jo + vsZ[i].dJ;}

		#pragma omp parallel for private(i,j,Jxyz,Jxzy,Jyxz,Jyzx,Jzxy,Jzyx) num_threads(num_threads)	  
		for (j=0;j<nc_plane;j++){ 
						Jxyz=0, Jxzy=0, Jyzx=0, Jyxz=0, Jzyx=0, Jzxy=0;
						if((ncx!=1)||(cube==1)){
								for (i=0;i<nsurx;i++){
										Jyxz+= vsX[i].J * emHyxz[j][i]; 						//Vi is included in interaction matrix emHyxz ...
										Jzxy+= vsX[i].J	* emHzxy[j][i];}}

						if((ncy!=1)||(cube==1)){
								for (i=0;i<nsury;i++){
										Jxyz+= vsY[i].J * emHxyz[j][i]; 
										Jzyx+= vsY[i].J * emHzyx[j][i];}}

						if((ncz!=1)||(cube==1)){
								for (i=0;i<nsurz;i++){
										Jxzy+= vsZ[i].J * emHxzy[j][i]; 
										Jyzx+= vsZ[i].J	* emHyzx[j][i];}}
 																								
						BJ_plane[j][0] = Jxyz - Jxzy;
						BJ_plane[j][1] = Jyzx - Jyxz; 
						BJ_plane[j][2] = Jzxy - Jzyx; 

						B_plane[j][0] = BJ_plane[j][0] + B_a[0];
						B_plane[j][1] = BJ_plane[j][1] + B_a[1]; 
						B_plane[j][2] = BJ_plane[j][2] + B_a[2];}	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void class_cube::magnetic_moment(int n, double m[3])
{ double x[3]={0}, J0[3]={0}, J[3]={0};

	interJ0_cell(n,J0);																							//interpolated J0

	vec_add2(J0,vc[n].J, J);																			//J=Jo+J(Jo+dJ)
	vec_scale(J,0.50);																							//J/2.0		

	vec_rot(vc[n].rc,J,x);																				//x=(rxJ)
	vec_scale1(x,vc[n].V/2,m);																		//1/2 * V * x
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::magnetic_moment(int n, double m[3])
{ double x[3]={0};

vec_rot(vc[n].rc,vc[n].J,x);
vec_scale1(x,vc[n].V/2.0,m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::fast_matrix_Hx()
{ int j,k[3]={0};
	double R=0,x1,x2,x3;
	k[0]=10;k[1]=10;k[2]=1;

//	if((cx>=cy)&&(cx>=cz))	{R=cx;}
//	if((cy>=cx)&&(cy>=cz))	{R=cy;}
//	if((cz>=cx)&&(cz>=cy))	{R=cz;}

	R=2*R;
	cout << "The analytical formula radius 2R " << R << " m" << endl;

	cout << "cmHyxz, cmHzxy ..." << endl;
	#pragma omp parallel for private(j,x1,x2,x3) num_threads(num_threads)
	for(j=0;j<nc;j++)
			{
//					x1=fabs(vc[j].rc[0]-vsX[0].rc[0]);
//					x2=fabs(vc[j].rc[1]-vsX[0].rc[1]);
//					x3=fabs(vc[j].rc[2]-vsX[0].rc[2]);

//					if((x1<R)&&(x2<R)&&(x3<R)){ //formula
//																			/* 4 cmHyxz*/	cmHyxz[j] = num_b(j, 0, k,0,2)/*magfvv(vc[j].rc, vsX[0].rc, vsX[0].size, vc[j].size, 2)*/;
//																			/* 5 cmHzxy*/	cmHzxy[j] = num_b(j, 0, k,0,1)/*magfvv(vc[j].rc, vsX[0].rc, vsX[0].size, vc[j].size, 1)*/;}
//														else
{		//approximation
																		/* 4 cmHyxz*/	cmHyxz[j] = Bvv(vc[j].rc, vsX[0].rc, vsX[0].Vi,2);
																		/* 5 cmHzxy*/	cmHzxy[j] = Bvv(vc[j].rc, vsX[0].rc, vsX[0].Vi,1);}}
									
								
	cout << "cmHxyz, cmHzyx ..." << endl;
	#pragma omp parallel for private(j,x1,x2,x3) num_threads(num_threads)
	for(j=0;j<nc;j++)
			{
//					x1=fabs(vc[j].rc[0]-vsY[0].rc[0]);
//					x2=fabs(vc[j].rc[1]-vsY[0].rc[1]);
//					x3=fabs(vc[j].rc[2]-vsY[0].rc[2]);

//					if((x1<R)&&(x2<R)&&(x3<R)){ //formula
//																			/* 1 cmHxyz*/	cmHxyz[j] = num_b(j, 0, k,1,2)/*magfvv(vc[j].rc, vsY[0].rc, vsY[0].size, vc[j].size, 2)*/;
//																			/* 6 cmHzyx*/	cmHzyx[j] = num_b(j, 0, k,1,0)/*magfvv(vc[j].rc, vsY[0].rc, vsY[0].size, vc[j].size, 0)*/;}
//														else
{		//approximation
																		/* 1 cmHxyz*/	cmHxyz[j]	= Bvv(vc[j].rc, vsY[0].rc, vsY[0].Vi,2);
																		/* 6 cmHzyx*/	cmHzyx[j] = Bvv(vc[j].rc, vsY[0].rc, vsY[0].Vi,0);}}												
		
	cout << "cmHxyz, cmHzyx ..." << endl;
	#pragma omp parallel for private(j,x1,x2,x3) num_threads(num_threads)
	for(j=0;j<nc;j++)
			{
//					x1=fabs(vc[j].rc[0]-vsZ[0].rc[0]);
//					x2=fabs(vc[j].rc[1]-vsZ[0].rc[1]);
//					x3=fabs(vc[j].rc[2]-vsZ[0].rc[2]);

//					if((x1<R)&&(x2<R)&&(x3<R)){ //formula
//																			/* 2 cmHxzy*/	cmHxzy[j] = num_b(j, 0, k,2,1)/*magfvv(vc[j].rc, vsZ[0].rc, vsZ[0].size, vc[j].size, 1)*/;
//																			/* 3 cmHyzx*/	cmHyzx[j] = num_b(j, 0, k,2,0)/*magfvv(vc[j].rc, vsZ[0].rc, vsZ[0].size, vc[j].size, 0)*/;}
//														else
{		//approximation
																		/* 2 cmHxzy*/	cmHxzy[j]	= Bvv(vc[j].rc, vsZ[0].rc, vsZ[0].Vi,1);
																		/* 3 cmHyzx*/	cmHyzx[j] = Bvv(vc[j].rc, vsZ[0].rc, vsZ[0].Vi,0);}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::full_matrix_Hx()
{ int i,j,k[3]={0};
	double a,R=0,x1,x2,x3;

	k[0]=5;k[1]=5;k[2]=1;

//	if((vc[0].size[0]>=vc[0].size[1])&&(vc[0].size[0]>=vc[0].size[2]))	{R=vc[0].size[0]+vc[0].size[0]/300;}
//	if((vc[0].size[1]>=vc[0].size[0])&&(vc[0].size[1]>=vc[0].size[2]))	{R=vc[0].size[1]+vc[0].size[1]/300;}
//	if((vc[0].size[2]>=vc[0].size[0])&&(vc[0].size[2]>=vc[0].size[1]))	{R=vc[0].size[2]+vc[0].size[2]/300;}

	R=2*R;
	cout << "The analytical formula radius 2R " << R << " m" << endl;

	cout << "mHyxz, mHzxy " << endl;
	#pragma omp parallel for private(i,j,x1,x2,x3) num_threads(num_threads)
	for(j=0;j<nc;j++)
			{
				for(i=0;i<nsurx;i++){

//					x1=fabs(vc[j].rc[0]-vsX[i].rc[0]);
//					x2=fabs(vc[j].rc[1]-vsX[i].rc[1]);
//					x3=fabs(vc[j].rc[2]-vsX[i].rc[2]);

//					if((x1<R)&&(x2<R)&&(x3<R)){	//formula
//																			/* 4 mHyxz*/ mHyxz[j][i] =/*num_b(j, i, k,0,2)*/magfvv(vc[j].rc, vsX[i].rc, vsX[i].size, vc[j].size, 2);
//																			/* 5 mHzxy*/ mHzxy[j][i] =/*num_b(j, i, k,0,1)*/magfvv(vc[j].rc, vsX[i].rc, vsX[i].size, vc[j].size, 1);}
//															else
{	//approximation
																		/* 4 mHyxz*/ mHyxz[j][i] = Bvv(vc[j].rc, vsX[i].rc, vsX[i].Vi,2);
																		/* 5 mHzxy*/ mHzxy[j][i] = Bvv(vc[j].rc, vsX[i].rc, vsX[i].Vi,1);}}}

	cout << "mHxyz, mHzyx " << endl;
	#pragma omp parallel for private(i,j,x1,x2,x3) num_threads(num_threads)
	for(j=0;j<nc;j++)
			{
				for(i=0;i<nsury;i++){		

//					x1=fabs(vc[j].rc[0]-vsY[i].rc[0]);
//					x2=fabs(vc[j].rc[1]-vsY[i].rc[1]);
//					x3=fabs(vc[j].rc[2]-vsY[i].rc[2]);
			
//					if((x1<R)&&(x2<R)&&(x3<R)){	//formula
//																			/* 1 mHxyz*/ mHxyz[j][i] = /*num_b(j, i, k,1,2)*/magfvv(vc[j].rc, vsY[i].rc, vsY[i].size, vc[j].size, 2);
//																			/* 6 mHzyx*/ mHzyx[j][i] = /*num_b(j, i, k,1,0)*/magfvv(vc[j].rc, vsY[i].rc, vsY[i].size, vc[j].size, 0);}
//															else
{	//approximation
																		/* 1 mHxyz*/ mHxyz[j][i] = Bvv(vc[j].rc, vsY[i].rc, vsY[i].Vi,2);
																		/* 6 mHzyx*/ mHzyx[j][i] = Bvv(vc[j].rc, vsY[i].rc, vsY[i].Vi,0);}}}

	if(ncz!=1){
	cout << "mHxzy, mHyzx " <<endl;
	#pragma omp parallel for private(i,j,x1,x2,x3) num_threads(num_threads)
	for(j=0;j<nc;j++)
			{
				for(i=0;i<nsurz;i++){	

//					x1=fabs(vc[j].rc[0]-vsZ[i].rc[0]);
//					x2=fabs(vc[j].rc[1]-vsZ[i].rc[1]);
//					x3=fabs(vc[j].rc[2]-vsZ[i].rc[2]);

//					if((x1<R)&&(x2<R)&&(x3<R)){	//formula
//																			/* 2 mHxzy*/ mHxzy[j][i] = /*num_b(j, i, k,2,1)*/magfvv(vc[j].rc, vsZ[i].rc, vsZ[i].size, vc[j].size, 1)*/;	
//																			/* 3 mHyzx*/ mHyzx[j][i] = /*num_b(j, i, k,2,0)*/magfvv(vc[j].rc, vsZ[i].rc, vsZ[i].size, vc[j].size, 0)*/;}
//															else
{	//approximation
																		/* 2 mHxzy*/ mHxzy[j][i] = Bvv(vc[j].rc, vsZ[i].rc, vsZ[i].Vi,1);
																		/* 3 mHyzx*/ mHyzx[j][i] = Bvv(vc[j].rc, vsZ[i].rc, vsZ[i].Vi,0);}}}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::matrix_eHx()
{ int i,j;
	double RC[3]={0};

	#pragma omp parallel for private(i,j) num_threads(num_threads)
	for(i=0;i<nc_plane;i++){
				for(j=0;j<nsurx;j++){
								/* 4 mHyxz */	emHyxz[i][j] = Bvv(RC_plane[i], vsX[j].rc, vsX[j].Vi,2);
								/* 5 mHzxy */	emHzxy[i][j] = Bvv(RC_plane[i], vsX[j].rc, vsX[j].Vi,1);}

				for(j=0;j<nsury;j++){
								/* 1 mHxyz */	emHxyz[i][j] = Bvv(RC_plane[i], vsY[j].rc, vsY[j].Vi,2);
								/* 6 mHzyx */ emHzyx[i][j] = Bvv(RC_plane[i], vsY[j].rc, vsY[j].Vi,0);}

				for(j=0;j<nsurz;j++){
								/* 2 mHxzy */ emHxzy[i][j] = Bvv(RC_plane[i], vsZ[j].rc, vsZ[j].Vi,1);
								/* 3 mHyzx */ emHyzx[i][j] = Bvv(RC_plane[i], vsZ[j].rc, vsZ[j].Vi,0);}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::read_cmHxyz(int i, int j)
{int jx=0,jy=0,jz=0,ix=0,iy=0,iz=0,kx=0,ky=0,kz=0;
	
	if(full_matrix==1){return mHxyz[i][j];}

  ix=vc[i].adr[0];				//		i=n/nsury;  i row cell received
  iy=vc[i].adr[1];
  iz=vc[i].adr[2];    

	jx=vsY[j].adr[0];				//		j=n%nsury;  j colum surface X created
	jy=vsY[j].adr[1];
	jz=vsY[j].adr[2];

	if(ix-jx>0) {kx=ix-jx;}else {kx=jx-ix;}
	if((2*(iy-jy)+1)>0) {ky=abs(iy-jy);}else {ky=abs(jy-iy-1);}
	if(iz-jz>0) {kz=iz-jz;}else {kz=jz-iz;}

	if(jz-iz==0) {return 0;}
	if(jz-iz<0) {return cmHxyz[icell(kx,ky,kz)];}
		 else{
				if(jz-iz>=0){return -cmHxyz[icell(kx,ky,kz)];}}						
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::read_cmHxzy(int i, int j)
{int jx=0,jy=0,jz=0,ix=0,iy=0,iz=0,kx=0,ky=0,kz=0;

	if(full_matrix==1){return mHxzy[i][j];}

  ix=vc[i].adr[0];						//		i=n/nsurz;  i row cell received
  iy=vc[i].adr[1];
  iz=vc[i].adr[2];    

	jx=vsZ[j].adr[0];				//		j=n%nsurz;  j colum surface X created
	jy=vsZ[j].adr[1];
	jz=vsZ[j].adr[2];

	if(ix-jx>0) {kx=ix-jx;}else {kx=jx-ix;}
	if(iy-jy>0) {ky=iy-jy;}else {ky=jy-iy;}
	if((2*(iz-jz)+1)>0) {kz=abs(iz-jz);}else {kz=abs(jz-iz-1);}

	if(jy-iy==0) {return 0;}
	if(jy-iy<0) {return cmHxzy[icell(kx,ky,kz)];}
		 else{
				if(jy-iy>0) {return -cmHxzy[icell(kx,ky,kz)];}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::read_cmHyzx(int i, int j)
{int jx=0,jy=0,jz=0,ix=0,iy=0,iz=0,kx=0,ky=0,kz=0;

	if(full_matrix==1){return mHyzx[i][j];}

  ix=vc[i].adr[0];						//		i=n/nsurz;  i row cell received
  iy=vc[i].adr[1];
  iz=vc[i].adr[2]; 

	jx=vsZ[j].adr[0];				//		j=n%nsury;  j colum surface X created
	jy=vsZ[j].adr[1];
	jz=vsZ[j].adr[2];

	if(ix-jx>0) {kx=ix-jx;}else {kx=jx-ix;}
	if(iy-jy>0) {ky=iy-jy;}else {ky=jy-iy;}
	if((2*(iz-jz)+1)>0) {kz=abs(iz-jz);}else {kz=abs(jz-iz-1);}

	if(jx-ix==0) {return 0;}
	if(jx-ix<0) {return cmHyzx[icell(kx,ky,kz)];}
		 else{
					if(jx-ix>0) {return -cmHyzx[icell(kx,ky,kz)];}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::read_cmHyxz(int i, int j)
{int jx=0,jy=0,jz=0,ix=0,iy=0,iz=0,kx=0,ky=0,kz=0;

	if(full_matrix==1){return mHyxz[i][j];}

  ix=vc[i].adr[0];						//		i=n/nsurx;  i row cell received
  iy=vc[i].adr[1];
  iz=vc[i].adr[2]; 

	jx=vsX[j].adr[0];				//		j=n%nsurx;  j colum surface X created
	jy=vsX[j].adr[1];
	jz=vsX[j].adr[2];



	if((2*(ix-jx)+1)>0) {kx=abs(ix-jx);}else {kx=abs(jx-ix-1);}
	if(iy-jy>0) {ky=iy-jy;}else {ky=jy-iy;}
	if(iz-jz>0) {kz=iz-jz;}else {kz=jz-iz;}

	if(jz-iz==0) {return 0;}
	if(jz-iz<0) {return cmHyxz[icell(kx,ky,kz)];}
		 else{
					if(jz-iz>0) {return -cmHyxz[icell(kx,ky,kz)];}}			
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::read_cmHzxy(int i, int j)
{int jx=0,jy=0,jz=0,ix=0,iy=0,iz=0,kx=0,ky=0,kz=0;

	if(full_matrix==1){return mHzxy[i][j];}

  ix=vc[i].adr[0];						//		i=n/nsurx;  i row cell received
  iy=vc[i].adr[1];
  iz=vc[i].adr[2]; 

	jx=vsX[j].adr[0];				//		j=n%nsurx;  j colum surface X created
	jy=vsX[j].adr[1];
	jz=vsX[j].adr[2];

	if ((2*(ix-jx)+1)>0) {kx=abs(ix-jx);}else {kx=abs(jx-ix-1);}
	if (iy-jy>0) {ky=iy-jy;}else {ky=jy-iy;}
	if (iz-jz>0) {kz=iz-jz;}else {kz=jz-iz;}

	if(jy-iy==0) {return 0;}
	if(jy-iy<0) {return cmHzxy[icell(kx,ky,kz)];}
		 else{
					if(jy-iy>0) {return -cmHzxy[icell(kx,ky,kz)];}}		
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::read_cmHzyx(int i, int j)
{ int jx=0,jy=0,jz=0,ix=0,iy=0,iz=0,kx=0,ky=0,kz=0;

	if(full_matrix==1){return mHzyx[i][j];}

  ix=vc[i].adr[0];						//		i=n/nsury;  i row cell received
  iy=vc[i].adr[1];
  iz=vc[i].adr[2]; 

	jx=vsY[j].adr[0];				//		j=n%nsury;  j colum surface X created
	jy=vsY[j].adr[1];
	jz=vsY[j].adr[2];

	if(ix-jx>0) {kx=ix-jx;}else {kx=jx-ix;}
	if((2*(iy-jy)+1)>0) {ky=abs(iy-jy);}else {ky=abs(jy-iy-1);}
	if(iz-jz>0) {kz=iz-jz;}else {kz=jz-iz;}

	if(jx-ix==0) {return 0;}
	if(jx-ix<0) {return cmHzyx[icell(kx,ky,kz)];}
		 else{
					if(jx-ix>0) {return -cmHzyx[icell(kx,ky,kz)];}}						
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::recal_theta()
{ int i;
				
	for(i=0;i<nc;i++){vc[i].theta = cell_theta(i);}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::cell_theta(int i)
{double theta, B;

	theta = atan2(vc[i].B[0],vc[i].B[2]) * 180 /pi;

	if(theta<0){theta = theta + 360;}

	return theta;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
								
double class_cube::num_b(int m, int l, int *na, int s, int c)			
{int i[3]={0},j[3]={0},A[3]={1,1,1},B;
 double x[3]={0}, x0[3]={0}, d[3]={0}, Bx=0, h, a[3]={0}, rc[3]={0};

	if(s==0){
	vec_cpy(a, vsX[l].size);														//size of surface a(x,y,z)
	vec_cpy(rc,vsX[l].rc);}															//position vector of surface l

	if(s==1){
	vec_cpy(a, vsY[l].size);														//size of surface a(x,y,z)
	vec_cpy(rc,vsY[l].rc);}				 											//position vector of surface l

	if(s==2){
	vec_cpy(a, vsZ[l].size);														//size of surface a(x,y,z)
	vec_cpy(rc,vsZ[l].rc);}															//position vector of surface l

	vec_div(d,a,na);																		//size of sub-element, na of sub-elements

	oper2(x0, rc, a, d, s);															//zero position of sub-element

	A[s]=2;

	for(i[0]=0;i[0]<na[0]*A[0];i[0]++){
			for(i[1]=0;i[1]<na[1]*A[1];i[1]++){
					for(i[2]=0;i[2]<na[2]*A[2];i[2]++){
								oper1(x, x0, i, d);																					//x0 + i*d change of sub-element position
								h=(-fabs(x[s] - rc[s]) + a[s])/a[s]; if(h<0){h=0;}					//calculates h function
								Bx+=h * magfvv(vc[m].rc, x, d, vc[m].size, c);}}}	
	return Bx;		
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Current density ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::JcB(double B[3], double Jc, double theta)
{double B0;


	if(rel==4){		
		B0=1e-3;

		if(vec_mod(B)<=B0)	{Jc=Jo - (Jo-Jo/3.0)*(vec_mod(B)/B0);}
		else{Jc=Jo/3;}
		return Jc;}

	if(rel==3) {return bilinear_inter_JcB(vec_mod(B), theta);}
	if(rel==2) {return (Jo/(pow(1.0 + (vec_mod(B)/Bo), m)));}
	return Jc;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::initial_Jc(double& min)
{int i,m;

	min = Jo;
	for(i=0;i<nc;i++){					
				if(vc[i].Jc!=Jol){ 
												{vc[i].Jc = JcB(vc[i].B, vc[i].Jc, vc[i].theta);}
												if(min>vc[i].Jc) 	min=vc[i].Jc;}}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::critical_Jc()
{int i;

	if(shape!=4){
		#pragma omp parallel for private(i) num_threads(num_threads)
		for (i=0;i<nc;i++) 	{vc[i].Jc = JcB(vc[i].B, vc[i].Jc, vc[i].theta);}}
			else{
						#pragma omp parallel for private(i) num_threads(num_threads)
						for (i=0;i<nc;i++) 	{if(vc[i].Jc!=Jol)vc[i].Jc = JcB(vc[i].B, vc[i].Jc, vc[i].theta);}}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::symm_dJ()
{int i,j,k,m;
double a;

	#pragma omp parallel for private(i,j,k,m) num_threads(num_threads)
	for(k=0;k<=ncz/2;k++){		
		for(j=0;j<=ncy/2;j++){
			for(i=0;i<=ncx/2;i++){

					if((ncx!=1)||(cube==1)){

								///// first quadrant/////
								m=ijktonsurX(i,j,k);
								if(j==ncy/2){vsX[m].dJ=0;}

								///// second quadrant /////
								vsX[vsX[m].adrss[0]].dJ = vsX[m].dJ;			
								///// third quadrant /////
								vsX[vsX[m].adrss[1]].dJ = (-1.0) * vsX[m].dJ;		
								///// fourth quadrant /////
								vsX[vsX[m].adrss[2]].dJ = (-1.0) * vsX[m].dJ;		

								if((cube==1)&&(k!=ncz/2)){
										///// fifth octet /////
										vsX[vsX[m].adrss[3]].dJ = vsX[m].dJ;			
										///// sixth octet /////
										vsX[vsX[m].adrss[4]].dJ = vsX[m].dJ;					
										///// seventh octet /////
										vsX[vsX[m].adrss[5]].dJ = (-1.0) * vsX[m].dJ;	
										///// eight octet /////
										vsX[vsX[m].adrss[6]].dJ = (-1.0) * vsX[m].dJ;}}


			 		if((ncy!=1)||(cube==1)){
								///// first quadrant/////
								m=ijktonsurY(i,j,k);
								if(i==ncx/2){vsY[m].dJ=0;}

								///// second quadrant /////
								vsY[vsY[m].adrss[0]].dJ = (-1.0) * vsY[m].dJ;		
								///// third quadrant /////
								vsY[vsY[m].adrss[1]].dJ = vsY[m].dJ;						
								///// fourth quadrant /////
								vsY[vsY[m].adrss[2]].dJ = (-1.0) * vsY[m].dJ;	

								if((cube==1)&&(k!=ncz/2)){
										///// fifth octet /////
										vsY[vsY[m].adrss[3]].dJ = vsY[m].dJ;				
										///// sixth octet /////
										vsY[vsY[m].adrss[4]].dJ = (-1.0) * vsY[m].dJ;}
										///// seventh octet /////
										vsY[vsY[m].adrss[5]].dJ = vsY[m].dJ;					
										///// eight octet /////
										vsY[vsY[m].adrss[6]].dJ = (-1.0) * vsY[m].dJ;}

							if((ncz!=1)||(cube==1)){
										///// first quadrant/////
										m=ijktonsurZ(i,j,k);
										if((i==ncx/2)||(j==ncy/2)){vsZ[m].dJ=0;}

										///// third quadrant /////
										vsZ[vsZ[m].adrss[0]].dJ = (-1.0) * vsZ[m].dJ;
										///// fifth octet /////
										vsZ[vsZ[m].adrss[1]].dJ = (-1.0) * vsZ[m].dJ;	
										///// seventh octet /////
										vsZ[vsZ[m].adrss[2]].dJ = vsZ[m].dJ;					
										///// second quadrant /////
										vsZ[vsZ[m].adrss[3]].dJ = (-1.0) * vsZ[m].dJ;	
										///// fourth quadrant /////
										vsZ[vsZ[m].adrss[4]].dJ = vsZ[m].dJ;					
										///// sixth octet /////
										vsZ[vsZ[m].adrss[5]].dJ = vsZ[m].dJ;				
										///// eight octet /////
										vsZ[vsZ[m].adrss[6]].dJ = (-1.0) * vsZ[m].dJ;}}}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::load_Jcbdata()
{int i,j,n;
 double a[5]={0},b=0;

	string str;
	fstream f;

	f.open("data/JcBload.txt",ios::in|ios::app);
	getline(f,str);			

	for(i=0;i<measured_fields;i++){
		for(n=0;n<measured_points;n++){
					for(j=0;j<5;j++) {f >> a[j];}
																		 data[i][n][0]  = n;				//s
																		 data[i][n][1]  = a[1];			//B
																		 data[i][n][2]  = a[2];			//theta
																		 data[i][n][3]  = a[3];			//Jc
														if(nB==1)data[i][n][4]  = a[4];}}		//n

	f.close();

	if(it==1){cout <<"Loaded measured data" << endl;
	cout <<"i  B[T]  theta Jc[A/m2]  n" << endl;
	for(i=0;i<measured_fields;i++){cout << i << endl;
	for (n=0;n<measured_points;n++) {cout << data[i][n][0] << " " << data[i][n][1] << " " << data[i][n][2] << " " << data[i][n][3] << " " << data[i][n][4] << endl;}cout << endl;}}
 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Dissipation factor ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::xdiss_UJ(double J[3], double Jc, double B[3], double Ns, int mm)
{
	double Jcpe1,Jcpa1;
	double B0;

	if (rel<4)	{return (1.0/(Ns+1.0)) * (Ec*Jc) * (pow ( (vec_mod(J) / Jc) , (Ns+1.0)));}

///// constant Jcpe, sharp curve Jcpa
	if (rel==4)	{
						B0=1e-3;
						Jcpa1=Jcpa;

						if(vec_mod(B)<=B0){Jcpa1=Jcpe - (Jcpe-Jcpa)*(vec_mod(B)/B0);}

						vc[mm].Jcpar = Jcpa1;

						if(vec_mod(B)<1e-10) return Uo*pow( vec_mod(J)/Jcpa1 , 2*mo );
						else return Uo*pow( (pow(vec_dot(J, B) / Jcpa1,2.0) + pow(vec_rot1(J, B) / Jcpe,2.0)) / pow(vec_mod(B),2.0), mo);}

/*
///// constant Jcpa, sharp curve Jcpe
	if (rel==4)	
				{
							B0=1e-3;
							Jcpe1=Jcpe;

							if(vec_mod(B)<=B0)	
								{
									Jcpe1=Jcpa - (Jcpa-Jcpe)*(vec_mod(B)/B0);
								}
	
						if(vec_mod(B)<1e-10) return Uo*pow( vec_mod(J)/Jcpa , 2*mo );
						else return Uo*pow( (pow(vec_dot(J, B) / Jcpa,2.0) + pow(vec_rot1(J, B) / Jcpe1,2.0)) / pow(vec_mod(B),2.0), mo);				
				}
*/
/*
///// Kim model 
	if (rel==2)	
				{
						if (vec_mod(B)==0)
								{
									Jcpe1=Jcpe;
									Jcpa1=Jcpe;
								}
								else
										{
//											Jcpe1 = Jcpe/(pow(1.0 + (sqrt(B[2]*B[2])/Bo), m));												//Jo=1.3e10 A/m2, Bo=20 mT, m=0.5,
//											Jcpa1 = Jcpe/(pow(1.0 + (sqrt(B[2]*B[2])/(9.0*Bo)), m));									//Jo=1.3e10 A/m2, Bo=9*20 mT, m=0.5,
												Jcpe1 = Jcpe/(pow(1.0 + (vec_mod(B)/Bo), m));															//Jo=1.3e10 A/m2, Bo=20 mT, m=0.5,
												Jcpa1 = Jcpe/(pow(1.0 + (vec_mod(B)/(9.0*Bo)), m));												//Jo=1.3e10 A/m2, Bo=9*20 mT, m=0.5,
										}

						if (vec_mod(B)==0) {return Uo*pow( (pow(vec_mod(J) / Jcpe1, 2.0)),mo);}
						return Uo*pow( (pow(vec_dot(J, B) / Jcpa1,2.0) + pow(vec_rot1(J, B) / Jcpe1,2.0)) / pow(vec_mod(B),2.0), mo);				
				}
*/
/*
///// constant Jcpa,Jcpe
	if (rel==4)	
				{
						Jcpe1=Jcpe;
						Jcpa1=Jcpa;

						if (vec_mod(B)==0) {return Uo*pow( (pow(vec_mod(J) / Jcpe, 2.0)),mo);}
						return Uo*pow( (pow(vec_dot(J, B) / Jcpa,2.0) + pow(vec_rot1(J, B) / Jcpe,2.0)) / pow(vec_mod(B),2.0), mo);				
				}
*/
	if(rel==5){
					if(vec_mod(J)<Jc) {return 0;}
						else{return rhoR*( pow(vec_mod(J),2.0) - pow(Jc,2.0) )/2.0;}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::ES_EJ(double J[3], int ADR, double Ns, double Jcx, double B[3], double E[3])
{double fac=0, modJ=0, Jpa=0, Jpe=0, Jpev[3]={0}, epa[3]={0}, epe[3]={0}, mod_B=0, Jcpe1, B0;

	modJ = vec_mod(J);

	if(rel<4){
				if(modJ==0) {vec_cpy(E,zero);}
						else{
								fac = Ec * (pow ((modJ/Jcx), Ns)); if (fac>1e300) {fac=1e300;};

								E[0] = fac * (J[0]/modJ); 
								E[1] = fac * (J[1]/modJ);
								E[2] = fac * (J[2]/modJ);}}
	
	if(rel==4){
				if(modJ==0) {vec_cpy(E,zero);}
					else{	
							mod_B = vec_mod(B);

							B0=1e-3;
							Jcpe1=Jcpe;

							if(mod_B<=B0)	{Jcpe1=Jcpa - (Jcpa-Jcpe)*(vec_mod(B)/B0);}

							///////////////////////////////////////////////////////

							if(mod_B<1e-10){
									Jpa=0.0;
									Jpe=vec_mod(J);
									vec_cpy(epa,zero);}
									else{
									Jpa = fabs(vec_dot(J, B))/mod_B;
									Jpe = vec_rot1(J, B)/mod_B; 
									vec_divi(B, mod_B ,epa);}																	 //	epa = B / |B|

							Jpev[0] = J[0] - Jpa*epa[0];
							Jpev[1] = J[1] - Jpa*epa[1];
							Jpev[2] = J[2] - Jpa*epa[2];

							fac = 2*mo*Uo*pow( pow(Jpa/Jcpa,2.0)+pow(Jpe/Jcpe1,2.0) ,(mo-1));

							E[0] = fac * (Jpa * epa[0]/pow(Jcpa,2.0) + Jpev[0]/pow(Jcpe1,2.0));		
							E[1] = fac * (Jpa * epa[1]/pow(Jcpa,2.0) + Jpev[1]/pow(Jcpe1,2.0));
							E[2] = fac * (Jpa * epa[2]/pow(Jcpa,2.0) + Jpev[2]/pow(Jcpe1,2.0));}}

	if(rel==5){
				if(modJ<Jcx) {vec_cpy(E,zero);}
				if(modJ>=Jcx){
								E[0] = rhoR * J[0]; 
								E[1] = rhoR * J[1];
								E[2] = rhoR * J[2];}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::L3(int m, double J[3])
{
	return  vc[m].V*(xdiss_UJ(J, vc[m].Jc, vc[m].B, vc[m].N, m) - vc[m].U);	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::U(int m)
{double J[3]={0};

	interJ_cell(m, J);

	return xdiss_UJ(J, vc[m].Jc, vc[m].B, vc[m].N, m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::U_dTz(int m)
{int m1,m2,m3,m4;

	m1 = veZ[m].adrc[0];
	m2 = veZ[m].adrc[1];
	m3 = veZ[m].adrc[2];
	m4 = veZ[m].adrc[3];

	vc[m1].U = U(m1);
	vc[m2].U = U(m2);
	vc[m3].U = U(m3);
	vc[m4].U = U(m4);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::U_dTx(int m)
{int m1,m2,m3,m4;

	m1 = veX[m].adrc[0];
	m2 = veX[m].adrc[1];
	m3 = veX[m].adrc[2];
	m4 = veX[m].adrc[3];

	vc[m1].U = U(m1);
	vc[m2].U = U(m2);
	vc[m3].U = U(m3);
	vc[m4].U = U(m4);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::U_dTy(int m)
{int m1,m2,m3,m4;

	m1 = veY[m].adrc[0];
	m2 = veY[m].adrc[1];
	m3 = veY[m].adrc[2];
	m4 = veY[m].adrc[3];

	vc[m1].U = U(m1);
	vc[m2].U = U(m2);
	vc[m3].U = U(m3);
	vc[m4].U = U(m4);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::UeZ(int m, double Jx1, double Jx2, double Jy1, double Jy2)
{int m1,m2,m3,m4;
 double J1[3]={0}, J2[3]={0}, J3[3]={0}, J4[3]={0};

	m1 = veZ[m].adrc[0];
	m2 = veZ[m].adrc[1];
	m3 = veZ[m].adrc[2];
	m4 = veZ[m].adrc[3];

	interJ_cell(m1, J1);		J1[0]+=(1.0/2.0)*Jx1;		J1[1]+=(1.0/2.0)*Jy1;
	interJ_cell(m2, J2); 		J2[0]+=(1.0/2.0)*Jx2;		J2[1]+=(1.0/2.0)*Jy1;
	interJ_cell(m3, J3);		J3[0]+=(1.0/2.0)*Jx2;		J3[1]+=(1.0/2.0)*Jy2;
	interJ_cell(m4, J4);		J4[0]+=(1.0/2.0)*Jx1;		J4[1]+=(1.0/2.0)*Jy2;

	if(shape==5){
			if(veZ[m].adr[1]==0)   {return 2*(L3(m2, J2) + L3(m3, J3));}
			if(veZ[m].adr[1]==ncy) {return 2*(L3(m1, J1) + L3(m4, J4));}}

	return L3(m1, J1) + L3(m2, J2) + L3(m3, J3) + L3(m4, J4);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::UeX(int m, double Jy1, double Jy2, double Jz1, double Jz2)
{int m1,m2,m3,m4;
 double J1[3]={0}, J2[3]={0}, J3[3]={0}, J4[3]={0};

	m1 = veX[m].adrc[0];
	m2 = veX[m].adrc[1];
	m3 = veX[m].adrc[2];
	m4 = veX[m].adrc[3];

	interJ_cell(m1, J1);		J1[1]+=(1.0/2.0)*Jy1;		J1[2]+=(1.0/2.0)*Jz2;
	interJ_cell(m2, J2);		J2[1]+=(1.0/2.0)*Jy2;		J2[2]+=(1.0/2.0)*Jz2;
	interJ_cell(m3, J3);		J3[1]+=(1.0/2.0)*Jy2;		J3[2]+=(1.0/2.0)*Jz1;
	interJ_cell(m4, J4);		J4[1]+=(1.0/2.0)*Jy1;		J4[2]+=(1.0/2.0)*Jz1;

	if(shape==5){
			if(veX[m].adr[1]==0)   {return 2*(L3(m1, J1) + L3(m2, J2));}
			if(veX[m].adr[1]==ncy) {return 2*(L3(m3, J3) + L3(m4, J4));}}

	return	L3(m1, J1) + L3(m2, J2) + L3(m3, J3) + L3(m4, J4);	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::UeY(int m, double Jx1, double Jx2, double Jz1, double Jz2)
{int m1,m2,m3,m4;
 double J1[3]={0}, J2[3]={0}, J3[3]={0}, J4[3]={0};

	m1 = veY[m].adrc[0];
	m2 = veY[m].adrc[1];
	m3 = veY[m].adrc[2];
	m4 = veY[m].adrc[3];

	interJ_cell(m1, J1);		J1[0]+=(1.0/2.0)*Jx1;		J1[2]+=(1.0/2.0)*Jz1;
	interJ_cell(m2, J2); 		J2[0]+=(1.0/2.0)*Jx2;		J2[2]+=(1.0/2.0)*Jz1;
	interJ_cell(m3, J3);		J3[0]+=(1.0/2.0)*Jx2;		J3[2]+=(1.0/2.0)*Jz2;
	interJ_cell(m4, J4);		J4[0]+=(1.0/2.0)*Jx1;		J4[2]+=(1.0/2.0)*Jz2;

	return	L3(m1, J1) + L3(m2, J2) + L3(m3, J3) + L3(m4, J4);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::recal_U()
{int m;

	#pragma omp parallel for private(m) num_threads(num_threads)
	for (m=0;m<nc;m++)	{vc[m].U = U(m);}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::recal_nB()
{int m;

	if(shape!=4){
	#pragma omp parallel for private(m) num_threads(num_threads)
	for (m=0;m<nc;m++)	{vc[m].N = bilinear_inter_nB(vec_mod(vc[m].B),vc[m].theta);}}
	else{
			#pragma omp parallel for private(m) num_threads(num_threads)
			for (m=0;m<nc;m++)	{if(vc[m].N!=1)vc[m].N = bilinear_inter_nB(vec_mod(vc[m].B),vc[m].theta);}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////// T vector //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::Jx_Tdl(int m)
{ int y1,y2,z1,z2;

	y1 = vsX[m].adre[0];
	y2 = vsX[m].adre[1];
	z1 = vsX[m].adre[2];
	z2 = vsX[m].adre[3];

	return (veY[y1].dT * veY[y1].l  +  veZ[z1].dT * veZ[z1].l  -  veY[y2].dT * veY[y2].l  -  veZ[z2].dT * veZ[z2].l)/vsX[m].S;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::Jy_Tdl(int m)
{ int x1,x2,z1,z2;

	x1 = vsY[m].adre[0];
	x2 = vsY[m].adre[1];
	z1 = vsY[m].adre[2];
	z2 = vsY[m].adre[3];

	return (-veX[x1].dT * veX[x1].l  -  veZ[z1].dT * veZ[z1].l  +  veX[x2].dT * veX[x2].l  +  veZ[z2].dT * veZ[z2].l)/vsY[m].S;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::Jz_Tdl(int m)
{ int x1,x2,y1,y2;

	x1 = vsZ[m].adre[0];
	x2 = vsZ[m].adre[1];
	y1 = vsZ[m].adre[2];
	y2 = vsZ[m].adre[3];

return (veX[x1].dT * veX[x1].l  +  veY[y1].dT * veY[y1].l  -  veX[x2].dT * veX[x2].l  -  veY[y2].dT * veY[y2].l)/vsZ[m].S;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::cal_dJfromdT()
{int i,j,k,m;

	if(sym==1){
				if(ncx!=1){
						#pragma omp parallel for private(m) num_threads(num_threads)
						for(m=0;m<nsurx;m++)	{vsX[m].dJ = Jx_Tdl(m);}}

				if(ncy!=1){
						#pragma omp parallel for private(m) num_threads(num_threads)
						for(m=0;m<nsury;m++)	{vsY[m].dJ = Jy_Tdl(m);}}

				if(ncz!=1){
						#pragma omp parallel for private(m) num_threads(num_threads)
						for(m=0;m<nsurz;m++)	{vsZ[m].dJ = Jz_Tdl(m);}}}

	if(sym==2){ 
					#pragma omp parallel for private(i,j,k,m) num_threads(num_threads)
					for(k=0;k<=ncz/2;k++){		
						for(j=0;j<=ncy/2;j++){
							for(i=0;i<=ncx/2;i++){	
										m = ijktonsurX(i,j,k);
										vsX[m].dJ = Jx_Tdl(m);

										m = ijktonsurY(i,j,k);
										vsY[m].dJ = Jy_Tdl(m);

										m = ijktonsurZ(i,j,k);
									 	vsZ[m].dJ = Jz_Tdl(m);}}}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::symm_dT()
{int i,j,k,m;

	#pragma omp parallel for private(i,j,k,m) num_threads(num_threads)
	for(k=0;k<=ncz/2;k++){			
		for(j=0;j<=ncy/2;j++){
			for(i=0;i<=ncx/2;i++){

					if((ncz==1)||(cube==1)){
								///// first quadrant/////
								m=ijktonedgeZ(i,j,k); 	

								///// second quadrant /////
								veZ[veZ[m].adres[0]].dT = veZ[m].dT;		
								///// third quadrant /////
  							veZ[veZ[m].adres[1]].dT = veZ[m].dT;		
								///// fourth quadrant /////
								veZ[veZ[m].adres[2]].dT = veZ[m].dT;		
						
								if((cube==1)&&(k!=ncz/2)){
										///// fifth octet /////
										veZ[veZ[m].adres[3]].dT = veZ[m].dT;
										///// sixth octet /////
										veZ[veZ[m].adres[4]].dT = veZ[m].dT; 
										///// seventh octet /////
										veZ[veZ[m].adres[5]].dT = veZ[m].dT;	
										///// eight octet /////
										veZ[veZ[m].adres[6]].dT = veZ[m].dT;}}


					if((ncx==1)||(cube==1)){
								///// first octet /////
								m=ijktonedgeX(i,j,k); 
								if(i==ncx/2){veX[m].dT=0;}

								///// third octet /////
								veX[veX[m].adres[1]].dT = veX[m].dT;						
								///// fifth octet /////
								veX[veX[m].adres[3]].dT = (-1.0) * veX[m].dT;
								///// seventh octet /////
								veX[veX[m].adres[5]].dT = (-1.0) * veX[m].dT;

								if((cube==1)&&(i!=ncx/2)){
											///// second octet /////
											veX[veX[m].adres[0]].dT = (-1.0) * veX[m].dT;
											///// fourth octet /////
											veX[veX[m].adres[2]].dT = (-1.0) * veX[m].dT;
											///// sixth octet /////
											veX[veX[m].adres[4]].dT = veX[m].dT;						
											///// eight octet /////	
											veX[veX[m].adres[6]].dT = veX[m].dT;}}

					if((ncy==1)||(cube==1)){
								///// first octet /////
								m=ijktonedgeY(i,j,k); 
								if(j==ncy/2){veY[m].dT=0;}
	
								///// second octet /////
								veY[veY[m].adres[0]].dT = veY[m].dT;							
								///// fifth octet /////
								veY[veY[m].adres[3]].dT = (-1.0) * veY[m].dT;	
								///// sixth octet /////
								veY[veY[m].adres[4]].dT = (-1.0) * veY[m].dT;		

								if((cube==1)&&(j!=ncy/2)){
											///// third octet /////
											veY[veY[m].adres[1]].dT = (-1.0) * veY[m].dT;
											///// fourth octet /////
											veY[veY[m].adres[2]].dT = (-1.0) * veY[m].dT;	
											///// seventh octet /////
											veY[veY[m].adres[5]].dT = veY[m].dT;			
											///// eight octet /////
											veY[veY[m].adres[6]].dT = veY[m].dT;}}}}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////  Minimization  /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::sector_iter(int s)
{int mzp=0, mzn=0, mxp=0, mxn=0, myp=0, myn=0, m,i,j,k,a=1,b=0;
 double dLxp=0, dLxn=0, dLyp=0, dLyn=0, dLzp=0, dLzn=0, dLx=0, dLy=0, dLz=0, stepJ=0, L=0;

	if((shape==5)&&(cij[s]==0))			{a=0;}
	if((shape==5)&&(caj[s]==ncy-1)) {b=1;}

		do{
			L=0;
			dLzp=1e300, dLzn=1e300;
			dLxp=1e300, dLxn=1e300;
			dLyp=1e300, dLyn=1e300;

			if((ncz==1)||(cube==1)){						
						for(k=cik[s];k<=cak[s];k++){	
								for(j=cij[s]+a;j<=caj[s]+b;j++){		
										for(i=cii[s]+1;i<=cai[s];i++){
														m = ijktonedgeZ(i,j,k);

														dLz =	dTz(m, htpz);
														if (dLz<dLzp) {dLzp=dLz; mzp=m;}	

														dLz = dTz(m, htnz);
														if (dLz<dLzn) {dLzn=dLz; mzn=m;}}}}}

			if((ncx==1)||(cube==1)){
						for(k=cik[s]+1;k<=cak[s];k++){						
								for(j=cij[s]+a;j<=caj[s]+b;j++){		
										for(i=cii[s];i<=cai[s];i++){
														m = ijktonedgeX(i,j,k);

														dLx =	dTx(m, htpx);
														if (dLx<dLxp) {dLxp=dLx; mxp=m;}	

														dLx = dTx(m, htnx);
														if (dLx<dLxn) {dLxn=dLx; mxn=m;}}}}}

			if((ncy==1)||(cube==1)){
						for(k=cik[s]+1;k<=cak[s];k++){							
								for(j=cij[s];j<=caj[s];j++){		
										for(i=cii[s]+1;i<=cai[s];i++){
														m = ijktonedgeY(i,j,k);

														dLy =	dTy(m, htpy);
														if (dLy<dLyp) {dLyp=dLy; myp=m;}	

														dLy = dTy(m, htny);
														if (dLy<dLyn) {dLyn=dLy; myn=m;}}}}}

			if((ncz==1)||(cube==1)){				
					if ((dLzp<dLzn)&&(dLzp<0)&&(dLzp<dLxp)&&(dLzp<dLxn)&&(dLzp<dLyp)&&(dLzp<dLyn)) {update_dTz(mzp, htpz, s);L+=1;}
					if ((dLzn<dLzp)&&(dLzn<0)&&(dLzn<dLxp)&&(dLzn<dLxn)&&(dLzn<dLyp)&&(dLzn<dLyn)) {update_dTz(mzn, htnz, s);L+=1;}}

			if((ncy==1)||(cube==1)){	
					if ((dLyp<dLyn)&&(dLyp<0)&&(dLyp<dLxp)&&(dLyp<dLxn)&&(dLyp<dLzp)&&(dLyp<dLzn)) {update_dTy(myp, htpy, s);L+=1;}
					if ((dLyn<dLyp)&&(dLyn<0)&&(dLyn<dLxp)&&(dLyn<dLxn)&&(dLyn<dLzp)&&(dLyn<dLzn)) {update_dTy(myn, htny, s);L+=1;}}

			if((ncx==1)||(cube==1)){	
					if ((dLxp<dLxn)&&(dLxp<0)&&(dLxp<dLyp)&&(dLxp<dLyn)&&(dLxp<dLzp)&&(dLxp<dLzn)) {update_dTx(mxp, htpx, s);L+=1;}
					if ((dLxn<dLxp)&&(dLxn<0)&&(dLxn<dLyp)&&(dLxn<dLyn)&&(dLxn<dLzp)&&(dLxn<dLzn)) {update_dTx(mxn, htnx, s);L+=1;}}

		stepJ++;
		}
//		while(stepJ<1);
		while ((L==1)&&(stepJ<1e7));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////  Energy  ///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::dTz(int m, double ht)
{int x1,x2,y1,y2,i,j,k;
 double Jx1, Jx2, Jy1, Jy2, L11=0;

	x1=veZ[m].adrs[0];
	x2=veZ[m].adrs[1];
	y1=veZ[m].adrs[2];
	y2=veZ[m].adrs[3];

	Jx1   =      	   ht * veZ[m].l/vsX[x1].S;
	Jx2 	= (-1.0) * ht * veZ[m].l/vsX[x2].S;

	Jy1   =        	 ht * veZ[m].l/vsY[y1].S;
	Jy2   = (-1.0) * ht * veZ[m].l/vsY[y2].S;

	///////// boundary conditions for infinite strip ////////////////////
	if(shape==5){if(veZ[m].adr[1]==0){
															return   	 L1_y(y1, Jy1) + 2*L1_x(x2, Jx2) + L1_y(y2, Jy2)
																			 + L2_y(y1, Jy1) + 2*L2_x(x2, Jx2) + L2_y(y2, Jy2)

																			 + (1.0/dt) * Jx2 * Jx2 * vsX[x2].Vi * vsX[x2].Vi * read_cmAy(y2,y1)//not the best... it's ok as long sample is square
																			 + (1.0/dt) * Jy1 * Jy2 * vsY[y1].Vi * vsY[y2].Vi * read_cmAy(y2,y1)

																			 + UeZ(m, 0, Jx2, Jy1, Jy2);}

													 if(veZ[m].adr[1]==ncy){
															return     2*L1_x(x1, Jx1) + L1_y(y1, Jy1) + L1_y(y2, Jy2)
																			 + 2*L2_x(x1, Jx1) + L2_y(y1, Jy1) + L2_y(y2, Jy2)

																			 + (1.0/dt) * Jx1 * Jx1 * vsX[x1].Vi * vsX[x1].Vi * read_cmAy(y2,y1)//not the best... it's ok as long sample is square
																			 + (1.0/dt) * Jy1 * Jy2 * vsY[y1].Vi * vsY[y2].Vi * read_cmAy(y2,y1)

																			 + UeZ(m, Jx1, 0, Jy1, Jy2);}}

	return     L1_x(x1, Jx1) + L1_y(y1, Jy1) + L1_x(x2, Jx2) + L1_y(y2, Jy2)
					 + L2_x(x1, Jx1) + L2_y(y1, Jy1) + L2_x(x2, Jx2) + L2_y(y2, Jy2)

					 + (1.0/dt) * Jx1 * Jx2 * vsX[x1].Vi * vsX[x2].Vi * read_cmAx(x2,x1)
					 + (1.0/dt) * Jy1 * Jy2 * vsY[y1].Vi * vsY[y2].Vi * read_cmAy(y2,y1)

					 + UeZ(m, Jx1, Jx2, Jy1, Jy2);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::dTx(int m, double ht)
{int y1,y2,z1,z2;
 double Jz1, Jz2, Jy1, Jy2;

	y1=veX[m].adrs[0];
	y2=veX[m].adrs[1];
	z1=veX[m].adrs[2];
	z2=veX[m].adrs[3];

	Jy1		=				 	 ht * veX[m].l/vsY[y1].S;
	Jy2 	=	(-1.0) * ht * veX[m].l/vsY[y2].S;

	Jz1 	=	(-1.0) * ht * veX[m].l/vsZ[z1].S;
	Jz2 	= 				 ht * veX[m].l/vsZ[z2].S;

	return		 L1_z(z1, Jz1) + L1_y(y1, Jy1) + L1_z(z2, Jz2) + L1_y(y2, Jy2)	
					 + L2_z(z1, Jz1) + L2_y(y1, Jy1) + L2_z(z2, Jz2) + L2_y(y2, Jy2)	

					 + (1.0/dt) * Jy1 * Jy2 * vsY[y1].Vi * vsY[y2].Vi * read_cmAy(y2,y1)
					 + (1.0/dt) * Jz1 * Jz2 * vsZ[z1].Vi * vsZ[z2].Vi * read_cmAz(z2,z1)

					 + UeX(m, Jy1, Jy2, Jz1, Jz2);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::dTy(int m, double ht)
{int x1,x2,z1,z2;
 double Jx1, Jx2, Jz1, Jz2;

	x1=veY[m].adrs[0];
	x2=veY[m].adrs[1];
	z1=veY[m].adrs[2];
	z2=veY[m].adrs[3];

	Jx1 = (-1.0) * ht * veY[m].l/vsX[x1].S;
	Jx2 = 				 ht * veY[m].l/vsX[x2].S;

	Jz1 = (-1.0) * ht * veY[m].l/vsZ[z1].S;
	Jz2 = 				 ht * veY[m].l/vsZ[z2].S;

return 		 L1_x(x1, Jx1) + L1_z(z1, Jz1) + L1_x(x2, Jx2) + L1_z(z2, Jz2)
				 + L2_x(x1, Jx1) + L2_z(z1, Jz1) + L2_x(x2, Jx2) + L2_z(z2, Jz2)

				 + (1.0/dt) * Jx1 * Jx2 * vsX[x1].Vi * vsX[x2].Vi * read_cmAx(x2,x1)
				 + (1.0/dt) * Jz1 * Jz2 * vsZ[z1].Vi * vsZ[z2].Vi * read_cmAz(z2,z1) 

				 + UeY(m, Jx1, Jx2, Jz1, Jz2);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::L1_x(int m, double h)
{
	return vsX[m].Vi * (h/dt) * (	vsX[m].av_dA + vsX[m].Vi * (h/2.0) * read_cmAx(m,m));  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::L2_x(int m, double h)
{
	return vsX[m].Vi * h * vsX[m].dAa/dt;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::L1_y(int m, double h)
{					
	return vsY[m].Vi * (h/dt) * (	vsY[m].av_dA + vsY[m].Vi * (h/2.0) * read_cmAy(m,m)); 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::L2_y(int m, double h)
{							
	return vsY[m].Vi * h * (vsY[m].dAa + vsY[m].dAs)/dt;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::L1_z(int m, double h)
{								
	return vsZ[m].Vi * (h/dt) * (vsZ[m].av_dA + vsZ[m].Vi * (h/2.0) * read_cmAz(m,m)); 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::L2_z(int m, double h)
{
	return vsZ[m].Vi * h * vsZ[m].dAa/dt;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////// Update energy  /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::update_dTz(int m, double ht, int s)
{int n,i,j,k, x1, x2, y1, y2, y3, y4;
 double ax1, ax2, ay1, ay2, Jx1, Jx2, Jy1, Jy2;						

	veZ[m].dT+= ht;

	x1=veZ[m].adrs[0];
	x2=veZ[m].adrs[1];
	y1=veZ[m].adrs[2];
	y2=veZ[m].adrs[3];

	Jx1   =      	   ht * veZ[m].l/(vsX[x1].S);
	Jx2 	= (-1.0) * ht * veZ[m].l/(vsX[x2].S);

	Jy1   =        	 ht * veZ[m].l/(vsY[y1].S);
	Jy2   = (-1.0) * ht * veZ[m].l/(vsY[y2].S);

	ax1 = Jx1 * vsX[x1].Vi;
	ax2 = Jx2 * vsX[x2].Vi;
	ay1 = Jy1 * vsY[y1].Vi; 
	ay2 = Jy2 * vsY[y2].Vi; 

	for (k=cik[s];k<=cak[s];k++){							
			for (j=cij[s];j<=caj[s];j++){		
					for (i=cii[s];i<=cai[s]+1;i++){

									n =  ijktonsurX(i,j,k);
									vsX[n].av_dA+= ax1 * read_cmAx(x1,n) + ax2 * read_cmAx(x2,n);}}}						//recalculated average potential of surface X	
 
	for (k=cik[s];k<=cak[s];k++){							
			for (j=cij[s];j<=caj[s]+1;j++){		
					for (i=cii[s];i<=cai[s];i++){

									n =  ijktonsurY(i,j,k);
									vsY[n].av_dA+= ay1 * read_cmAy(y1,n) + ay2 * read_cmAy(y2,n);}}}						//recalculated average potential of surface Y

	vsX[x1].dJ+= Jx1;			
	vsX[x2].dJ+= Jx2;			
	vsY[y1].dJ+= Jy1;			
	vsY[y2].dJ+= Jy2;

	U_dTz(m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::update_dTx(int m, double ht, int s)
{int n,i,j,k, y1, y2, z1, z2, y3,y4;
 double az1, az2, ay1, ay2, Jz1, Jz2, Jy1, Jy2;

	veX[m].dT+= ht;

	y1=veX[m].adrs[0];
	y2=veX[m].adrs[1];
	z1=veX[m].adrs[2];
	z2=veX[m].adrs[3];

	Jy1 = 				 ht * veX[m].l/vsY[y1].S;
	Jy2 = (-1.0) * ht * veX[m].l/vsY[y2].S;

	Jz1 = (-1.0) * ht * veX[m].l/vsZ[z1].S;
	Jz2 = 				 ht * veX[m].l/vsZ[z2].S;

	////////////// Boundary strip /////////////////
	if(shape==5){if(veX[m].adr[1]==0)		{Jz1=0;}
							 if(veX[m].adr[1]==ncy) {Jz2=0;}}

	ay1 = Jy1 * vsY[y1].Vi; 
	ay2 = Jy2 * vsY[y2].Vi; 
	az1 = Jz1 * vsZ[z1].Vi;
	az2 = Jz2 * vsZ[z2].Vi;

	for (k=cik[s];k<=cak[s];k++){							
			for (j=cij[s];j<=caj[s]+1;j++){		
					for (i=cii[s];i<=cai[s];i++){

									n = ijktonsurY(i,j,k);
									vsY[n].av_dA+= ay1 * read_cmAy(y1,n) + ay2 * read_cmAy(y2,n);}}}				//recalculated average potential of surface Y
															
	for (k=cik[s];k<=cak[s]+1;k++){							
			for (j=cij[s];j<=caj[s];j++){		
					for (i=cii[s];i<=cai[s];i++){

									n = ijktonsurZ(i,j,k);
									vsZ[n].av_dA+= az1 * read_cmAz(z1,n) + az2 * read_cmAz(z2,n);}}}				//recalculated average potential of surface Z	

	vsY[y1].dJ+= Jy1;			
	vsY[y2].dJ+= Jy2;			
	vsZ[z1].dJ+= Jz1;			
	vsZ[z2].dJ+= Jz2;			

	U_dTx(m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::update_dTy(int m, double ht, int s)
{int n,i,j,k, x1, x2, z1, z2;
 double ax1, ax2, az1, az2, Jx1, Jx2, Jz1, Jz2;

	veY[m].dT+= ht;

	x1=veY[m].adrs[0];
	x2=veY[m].adrs[1];
	z1=veY[m].adrs[2];
	z2=veY[m].adrs[3];

	Jx1 = (-1.0) *	ht * veY[m].l/vsX[x1].S;
	Jx2 = 					ht * veY[m].l/vsX[x2].S;

	Jz1 = (-1.0) *	ht * veY[m].l/vsZ[z1].S;
	Jz2 =						ht * veY[m].l/vsZ[z2].S;

	ax1 = Jx1 * vsX[x1].Vi;
	ax2 = Jx2 * vsX[x2].Vi;
	az1 = Jz1 * vsZ[z1].Vi; 
	az2 = Jz2 * vsZ[z2].Vi; 

	for (k=cik[s];k<=cak[s];k++){							
			for (j=cij[s];j<=caj[s];j++){		
					for (i=cii[s];i<=cai[s]+1;i++){

									n =  ijktonsurX(i,j,k);
									vsX[n].av_dA+= ax1 * read_cmAx(x1,n) + ax2 * read_cmAx(x2,n);}}}				//recalculated average potential of surface X	


	for (k=cik[s];k<=cak[s]+1;k++){							
			for (j=cij[s];j<=caj[s];j++){		
					for (i=cii[s];i<=cai[s];i++){
									n =  ijktonsurZ(i,j,k);
									vsZ[n].av_dA+= az1 * read_cmAz(z1,n) + az2 * read_cmAz(z2,n);}}}				//recalculated average potential of surface Z	


	vsX[x1].dJ+= Jx1;			
	vsX[x2].dJ+= Jx2;		
	vsZ[z1].dJ+= Jz1;			
	vsZ[z2].dJ+= Jz2;			

	U_dTy(m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////  AC loss  //////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::local_loss(int m)
{double modJ,J0[3]={0},J[3]={0};

	interJ0_cell(m,J0);

	vec_add2(J0,vc[m].J, J);																				//Jo+J(Jo+dJ)
	vec_scale(J,0.50);																								//J/2.0		

	modJ = vec_mod(J);
	 
	if(modJ==0) {vec_cpy(vc[m].uv,zero);}													//set unit vector = 0						
			else{vec_divi(J,modJ, vc[m].uv);}													//unit vector=J/modJ

	ES_EJ(vc[m].J,m, vc[m].N, vc[m].Jc, vc[m].B, vc[m].E);

	vc[m].loss = vc[m].V * vec_dot(vc[m].E,J);//if(vc[m].N==Nl)cout<<vc[m].loss<<" "<<vc[m].E[0]<<" "<<vc[m].J[0]<<endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double class_cube::loss()
{int n,i,j,k,s;
 double J[3]={0}, m[3]={0}, m_add[3]={0}, a=1e-20, ff=0, Pm=0, XX=0, M[3]={0};

	ff=f;

	fstream f;
	f.open("data/loss.txt",ios::out|ios::app);
	if (it==1){
					f<<"#1  2  3		  4		    5		   6 	  7	   8	  9	  10  11  12  13  14	15    16     17 	 " <<endl;
					f<<"#n  t  Px[W]  t-dt/2  Ph[W]  Bax  Bay  Baz  mx  my  mz  Mx  My  Mz  M*ea  Ps[W]  Pl[W] " <<endl;
					f<<it-1<<"\t"<<tp<<"\t"<<Px<<"\t"<<t-0.5*dt<<"\t"<<Ph<<"\t"<<pBa[0]<<"\t"<<pBa[1]<<"\t"<<pBa[2]<<"\t"<<m[0]<<"\t"<<m[1]<<"\t"<<m[2]<<"\t"<<M[0]<<"\t"<<
						 M[1]<<"\t"<<M[2]<<"\t"<<XX<<"\t"<<Psup<<"\t"<<Plin<< endl;
					vec_cpy(mp,zero);}

	if ((it==11)&&((Btrape==1))){
					f<<"#1  2  3		  4		    5		   6 	  7	   8	  9	  10  11  12  13  14	15    16     17    18	 19	 20	 21 22 23" <<endl;
					f<<"#n  t  Px[W]  t-dt/2  Ph[W]  Bax  Bay  Baz  mx  my  mz  Mx  My  Mz  M*ea  Ps[W]  Pl[W] Rcx Rcy Rcz Bx By Bz" <<endl;
					f<<it-1<<"\t"<<tp<<"\t"<<Px<<"\t"<<t-1.5*dt<<"\t"<< 0 <<"\t"<<pBa[0]<<"\t"<<pBa[1]<<"\t"<<pBa[2]<<"\t"<<m[0]<<"\t"<<m[1]<<"\t"<<m[2]<<"\t"<<M[0]<<"\t"<<
						 M[1]<<"\t"<<M[2]<<"\t"<<XX<<"\t"<<Psup<<"\t"<<Plin<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<< endl;
					vec_cpy(mp,zero);}
	f.close();

	Pix=Px;
	Pilin=Plin;
	Pisup=Psup;

	Px=0; Ph=0; Plin=0; Psup=0;
	 
	for (n=0;n<nc;n++){
				Ph+=vc[n].loss;
				m[0]+=vc[n].m[0];
				m[1]+=vc[n].m[1];
				m[2]+=vc[n].m[2];}

	cout << " mx " << m[0] << " my " << m[1] << " mz " << m[2] << endl;

	if(fabs(m[0])<a) m[0]=0;
	if(fabs(m[1])<a) m[1]=0;
	if(fabs(m[2])<a) m[2]=0;

	vec_divi(m,V_tot,M);

	XX=vec_dot(M, ea);


	if(Btrape==1){s=iplanecell(ncx_plane/2, ncy_plane/2, 0);}	

	cout << "Magnetic field plane [position of the cell] " << ncx_plane/2 << " " << ncy_plane/2 << endl;

	f.open("data/loss.txt",ios::out|ios::app);

	if(((it*dt/T)>=0.25)&&((it*dt/T)<=1.25)){
															for(n=0;n<nc;n++){
																		Px+=vc[n].loss;	
																		if(vc[n].N==Nl) {Plin+=vc[n].loss;}
																		if(vc[n].N==N)  {Psup+=vc[n].loss;}}

															vec_add2(m, mp, m_add);					
															Pm=vec_dot(m_add, dBa);
															Qh+=-Pm/2.0;
															if((it*dt/T)==0.75){Qh1=Qh;}

															if((it*dt/T)>0.25){
																		Q+=  (-tp + t)*(Pix + Px)     * (1.0/2.0);
																		Ql+= (-tp + t)*(Pilin + Plin) * (1.0/2.0); 
																		Qs+= (-tp + t)*(Pisup + Psup) * (1.0/2.0);}}
																				else{for (n=0;n<nc;n++){
																										if(vc[n].N==Nl) {Plin+=vc[n].loss;}
																										if(vc[n].N==N)  {Psup+=vc[n].loss;}}}

	if(Btrape==1){
				f<< it <<"\t"<< setprecision(9) << t <<"\t" << Px << "\t" << t-0.5*dt << "\t" << Ph << "\t" << B_a[0] << "\t" << B_a[1] << "\t" << B_a[2] << "\t" << m[0] << "\t" << m[1] << "\t" << m[2] << "\t" << 
						M[0] << "\t" << M[1] << "\t" << M[2] << "\t" << XX <<"\t" << Psup << "\t" << Plin <<"\t"<< RC_plane[s][0] <<"\t"<< RC_plane[s][1] <<"\t"<< RC_plane[s][2] <<"\t"<< 
						B_plane[s][0] <<"\t"<< B_plane[s][1] <<"\t"<< B_plane[s][2] <<endl;}
			else{
				f<< it <<"\t"<< setprecision(9) << t <<"\t" << Px << "\t" << t-0.5*dt << "\t" << Ph << "\t" << B_a[0] << "\t" << B_a[1] << "\t" << B_a[2] << "\t" << m[0] << "\t" << m[1] << "\t" << m[2] << "\t" << 
						M[0] << "\t" << M[1] << "\t" << M[2] << "\t" << XX <<"\t" << Psup << "\t" << Plin <<"\t"<<endl;}

	f.close();
	vec_cpy(mp,m);	
}	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////// Is current  //////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::set_us()
{int i;

	#pragma omp parallel for private(i) num_threads(num_threads)
	for(i=0;i<nsury;i++){vsY[i].us=	1.0/(x*z);}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////  Other  ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_matrix()
{

	cout << "Cartezian matrix..." << endl;
	if(full_matrix==0)fast_matrix_cmA();
	if(full_matrix==1)full_matrix_mA();

	cout<<"Matrix Hx, Hy, Hz..." << endl;
	if(full_matrix==0)fast_matrix_Hx();
	if(full_matrix==1)full_matrix_Hx();

	if(Btrape==1){
			cout<<"Matrix eHx, eHy, eHz..." << endl;
			matrix_eHx();}

	if (nncx==0){cout << "Without linear material." << endl;}
	if (nncx!=0){cout << "With linear material." << endl;}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_temp()
{int m;

		if(ncx!=1){ 
					#pragma omp parallel for private(m) num_threads(num_threads)
					for(m=0;m<nsurx;m++)	{vsX[m].dJx = vsX[m].dJ;}}

		if(ncy!=1){
					#pragma omp parallel for private(m) num_threads(num_threads)
					for(m=0;m<nsury;m++)	{vsY[m].dJx = vsY[m].dJ;}}

		if(ncz!=1){
					#pragma omp parallel for private(m) num_threads(num_threads)
					for(m=0;m<nsurz;m++)	{vsZ[m].dJx = vsZ[m].dJ;}}

		if((ncx==1)||(cube==1)){
							#pragma omp parallel for private(m) num_threads(num_threads)
							for(m=0;m<nedgex;m++) {veX[m].dTx = veX[m].dT;}}

		if((ncy==1)||(cube==1)){
							#pragma omp parallel for private(m) num_threads(num_threads)
							for(m=0;m<nedgey;m++)	{veY[m].dTx = veY[m].dT;}}

		if((ncz==1)||(cube==1)){
					#pragma omp parallel for private(m) num_threads(num_threads)
					for(m=0;m<nedgez;m++)	{veZ[m].dTx = veZ[m].dT;}}
}	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::damping(double l)
{int m;

	if((ncz==1)||(cube==1)){
			#pragma omp parallel for private(m) num_threads(num_threads)
			for(m=0;m<nedgez;m++){veZ[m].dT = veZ[m].dTx + l * (veZ[m].dT - veZ[m].dTx);}}

	if((ncx==1)||(cube==1)){
			#pragma omp parallel for private(m) num_threads(num_threads)
			for (m=0;m<nedgex;m++){veX[m].dT = veX[m].dTx + l * (veX[m].dT - veX[m].dTx);}}

	if((ncy==1)||(cube==1)){
			#pragma omp parallel for private(m) num_threads(num_threads)
			for (m=0;m<nedgey;m++){veY[m].dT = veY[m].dTx + l * (veY[m].dT - veY[m].dTx);}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::dif(double& max, double l, int& sum)
{int m;
 double a=0,b=0;

	sum=0;
	max=0;

	if((ncz==1)||(cube==1)){
				for(m=0;m<nedgez;m++){
								b = veZ[m].dTx + l * (veZ[m].dT - veZ[m].dTx);
								a	= fabs(veZ[m].dTx - b);
								if (a>max) {max=a;}
								if(veZ[m].dT!=veZ[m].dTx) {sum++;}}}

	if((ncx==1)||(cube==1)){
					for(m=0;m<nedgex;m++){
								b = veX[m].dTx + l * (veX[m].dT - veX[m].dTx);
								a = fabs(veX[m].dTx - b);
								if (a>max) {max=a;}
								if(veX[m].dT!=veX[m].dTx) {sum++;}}}

	if((ncy==1)||(cube==1)){
					for(m=0;m<nedgey;m++){
								b = veY[m].dTx + l * (veY[m].dT - veY[m].dTx);
								a = fabs(veY[m].dTx - b);
								if (a>max) {max=a;}
								if(veY[m].dT!=veY[m].dTx) {sum++;}}}						
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::set_Ba(int iit)
{int itp,i;
 double a,tx,Ta,nx;
 it=iit;

	vec_cpy(pBa,B_a);

	cout << "time step " << it << endl;

	if(Bshape==0){t  = dt*it;	
								tp = dt*(it-1);
								Ba = Bamax * sin (2.0 * pi * f * t);}

	if(Bshape==1){ 				
					if(it==11){
											pBa[0]=0;
											pBa[1]=0;
											pBa[2]=Bamax;
											Ba=Bamax;
											time_variable();}

					t  = dt*it;	
					tp = dt*(it-1);
					itp=ns/4.0;
	
					if(t<=T){	////magnetize////
										////monopolar////
										if(t<=T/4.0) 	{Ba+= Bamax/itp;}																	//ramp up
										if((t>T/4.0)&&(t<=(T/2.0)*(1.0+1e-15))) {Ba-= Bamax/itp;}				//ramp down


										////bipolar////
		//								if(t<=T/4.0) 	{Ba+= Bamax/itp;}
		//								if((t>T/4)  &&(t<=(3*T/4)*(1.0+1e-15))) {Ba-= Bamax/itp;}
		//								if((t>3*T/4)&&(t<=T*(1.0+1e-15))) 		{Ba+= Bamax/itp;}
									}

					Ta=T*1.25;
					nx=ns*1.25;
					a=1.0+1e-5;


					if(it>nx){
									T1  = 1.0/f1;
									dt1 = T1/ns; 
									dt  = dt1;
									t   = Ta + (it-nx)*dt1;

									fi=fi1; 
									Bamax=Bamax1;
						
									Ba = Bamax * sin (2.0 * pi * f1 * (t-Ta));}}

	if(Bshape==2){
					a  = 1.0 + 1e-5;
					t  = dt*it;	
					tp = dt*(it-1);
					itp= ns/4.0;

					tx = (t-T)/T - int((t-T)/T);

					if(t<=T){tx = t/T - int(t/T);}
						else{tx = (t-T)/T - int((t-T)/T);}


					////bipolar////
					if(tx<=0.25*a) 								{Ba+= Bamax/itp;}
					if((tx>0.25*a)&&(tx<=0.75*a)) {Ba-= Bamax/itp;}
					if((tx>0.75*a)&&(tx<=1.0*a)) 	{Ba+= Bamax/itp;}}

	B_a[0] = Ba * sin(fi*pi/180.0)*cos(theta*pi/180.0);				if(fabs(B_a[0])<1e-12) {B_a[0]=0.0;}  //instant magnetic field B_a x component		
	B_a[1] = Ba * sin(fi*pi/180.0)*sin(theta*pi/180.0);				if(fabs(B_a[1])<1e-12) {B_a[1]=0.0;}  //instant magnetic field B_a y component
	B_a[2] = Ba * cos(fi*pi/180.0);											  		if(fabs(B_a[2])<1e-12) {B_a[2]=0.0;}  //instant magnetic field B_a z component

	vec_subs2(B_a,pBa,dBa);      
	cout << " Ba " << B_a[0]*1000 << " " << B_a[1]*1000 << " " << B_a[2]*1000 << "[mT]  dBa " << dBa[0]*1000 <<" "<< dBa[1]*1000 <<" "<< dBa[2]*1000 << " [mT]" << endl;	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::set_Is(int iit)
{int i;
  it=iit;

	pI=Is;

	if(Bamax==0){cout << "time step " << it << endl;}

	t  = dt*it;	
	tp = dt*(it-1);
	Is = Ismax * sin (2.0 * pi * f * t);
	
	dI=Is-pI;	
	cout << "Is " << Is << " [A] dI " << dI << " [A]" <<endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::set_Jc(int iit)
{int itp,i;
 double a,tx,Ta,nx;
 it=iit;

	if((rel==2)||(rel==3)){magnetic_field();}																			//magnetic field B= BA  

	#pragma omp parallel for private(i) num_threads(num_threads)	  
	for(i=0;i<nc;i++){vec_cpy(vc[i].pB,vc[i].B);}

	if(rel==3){load_Jcbdata();}																										//load Jc(B) data from measurement

	if(rel==1){minJc=Jo;}
				else{initial_Jc(minJc);}																								//crititcal current density according Kim model

	if(it!=1)recal_U();																														//recalculation of U in all edges X,Y,Z

	hp	 =          minJc; 
	hn	 = (-1.0) * minJc;

	if((ncz==1)||(cube==1)){
			if(x>y)
				{htpz	= hp * cy;																												// /10;/5;/2;????
				 htnz	= hn * cy;}																												// /10;/5;/2;????
				else
						{htpz = hp * cx;
						 htnz	= hn * cx;}
						 htpx=0;htpy=0;htnx=0;htny=0;}

	if((ncx==1)||(cube==1)){
			if(y>z)
				{htpx = hp * cz;
				 htnx = hn * cz;}
				else
						{htpx = hp * cy;
						 htnx = hn * cy;}}

	if((ncx==1)&&(cube==0)){
					htpx = hp * cy;
 				  htnx = hn * cy;
					htpy=0;htpz=0;htny=0;htnz=0;}


	if((ncy==1)||(cube==1)){
			if(x>z)
				{htpy = hp * cz;
				 htny	= hn * cz;}
				else
						{htpy	= hp * cx;
						 htny	= hn * cx;}}																									

	cout << "minJc " << minJc << endl;
	cout << "htpz " << htpz << endl;
	cout << "htnz " << htnz << endl << flush;
	cout << "htpx " << htpx << endl;
	cout << "htnx " << htnx << endl << flush;
	cout << "htpy " << htpy << endl;
	cout << "htny " << htny << endl << flush;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::time_variable() 
{int n,i,j,k;

	#pragma omp parallel for private(n) num_threads(num_threads)
	for(n=0;n<nsurx;n++){
          vsX[n].Jo+=vsX[n].dJ;
			  	vsX[n].dJ = 0;
			 		vsX[n].Aa0 = vsX[n].Aa0 + vsX[n].dAa;

					vsX[n].dAa = surX_dA(n);
					vsX[n].av_dA = 0;}

	#pragma omp parallel for private(n) num_threads(num_threads)
	for(n=0;n<nsury;n++){
          vsY[n].Jo+= vsY[n].dJ;
		    	vsY[n].dJ = 0;
			 		vsY[n].Aa0 = vsY[n].Aa0 + vsY[n].dAa;

					vsY[n].dAa = surY_dA(n);
					vsY[n].av_dA = 0;			

					if(Ismax!=0){
							vsY[n].dI=dI;
							vsY[n].Is=Is;
							vsY[n].dAs = vsY[n].dI * snAy[n];}}

	#pragma omp parallel for private(n) num_threads(num_threads)
	for(n=0; n<nsurz; n++){				
          vsZ[n].Jo+= vsZ[n].dJ;
	      	vsZ[n].dJ = 0;
			 		vsZ[n].Aa0 = vsZ[n].Aa0 + vsZ[n].dAa;

					vsZ[n].dAa = surZ_dA(n);
					vsZ[n].av_dA = 0;}

	if(Ismax!=0){
	#pragma omp parallel for private(n) num_threads(num_threads)
	for(n=0; n<nc; n++){				
          vc[n].Js[0] = 0;
          vc[n].Js[1] = interJs_cell(n, 1);
          vc[n].Js[2] = 0;}}

	#pragma omp parallel for private(n) num_threads(num_threads)
	for(n=0; n<nedgez; n++){
					veZ[n].T0+=veZ[n].dT;
					veZ[n].dT=0;}

	if(ncz!=1){
				#pragma omp parallel for private(n) num_threads(num_threads)
				for(n=0; n<nedgex; n++){
								veX[n].T0+=veX[n].dT;
								veX[n].dT=0;}

				#pragma omp parallel for private(n) num_threads(num_threads)
				for(n=0; n<nedgey; n++){
								veY[n].T0+=veY[n].dT;
								veY[n].dT=0;}}
}	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////  Save to file  /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::clear_files()
{
	fstream f;
	f.open("data/nodes.txt",ios::out);	f.clear();	f.close();
	f.open("data/edgesX.txt",ios::out);	f.clear();	f.close();
	f.open("data/edgesY.txt",ios::out);	f.clear();	f.close();
	f.open("data/edgesZ.txt",ios::out);	f.clear();	f.close();
	f.open("data/surfaces_X.txt",ios::out);	f.clear();	f.close();
	f.open("data/surfaces_Y.txt",ios::out);	f.clear();	f.close();
	f.open("data/surfaces_Z.txt",ios::out);	f.clear();	f.close();
	f.open("data/cells.txt",ios::out);	f.clear();	f.close();
	f.open("data/output_Jx.txt",ios::out);	f.clear();	f.close();
	f.open("data/output_Jy.txt",ios::out);	f.clear();	f.close();
	f.open("data/output_Jz.txt",ios::out);	f.clear();	f.close();
	f.open("data/loss.txt",ios::out);	f.clear();	f.close();
	f.open("data/output3Dx.txt",ios::out);	f.clear();	f.close();
	f.open("data/output3Dy.txt",ios::out);	f.clear();	f.close();
	f.open("data/output3Dz.txt",ios::out);	f.clear();	f.close();
	f.open("data/cross_sectionX.txt",ios::out);	f.clear();	f.close();
	f.open("data/cross_sectionY.txt",ios::out);	f.clear();	f.close();
	f.open("data/Tz.txt",ios::out);	f.clear();	f.close();
	f.open("data/outputAV.txt",ios::out);	f.clear();	f.close();
	f.open("data/output_plane.txt",ios::out);	f.clear();	f.close();
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_par()
{double ff=f;

	fstream f;
	f.open("data/par_plot.plt",ios::out);

	f << "x= " << x << endl;
	f << "y= " << y << endl;
	f << "z= " << z << endl;
	f << "ncx= " << ncx << endl; 
	f << "ncy= " << ncy << endl; 
	f << "ncz= " << ncz << endl; 
	f << "step= " << step << endl;
	f << "Bamax= " << Bamax << endl;
	f << "Ismax= " << Ismax << endl;
	f << "ns= " << ns << endl;
	f << "N= " << N << endl;
	f << "T= " << T << endl;
	f << "f= " << ff << endl;
	f << "dt= " << dt << endl; 
	f << "Jc= " << Jo << endl;
	f << "Jcpa= " << Jcpa << endl;
	f << "Jcpe= " << Jcpe << endl;
	f << "rel= " << rel << endl;

	f.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_3d()
{int i,j,k,m, n=0;
 double J[3]={0}, Jc1=0,Nfac=0, A[3]={0}, B[3]={0}, B_J[3]={0}, T[3]={0}, mm[3]={0}, P=0, e[3]={0},E[3]={0},U=0, Jcpar=0,I;

	fstream f;
	f.open("data/output3Dz.txt",ios::out|ios::app);
	if (it==1) {f <<"#1 2     3     4 		5 		 6 			 7  		8 		9  10 11 12  13  14  15  16   17   18 	19   20 	21 	 22 23 24 25 26 27 28 29 30 31 32 33 34 	35	36	37	 38"<< endl;
							f <<"#n rc[x] rc[y] rc[z] rec[r] rec[fi] rec[z] Jcpar Jx Jy Jz |J| Jc  U   N   A_ax A_ay A_az B_Jx B_Jy B_Jz Bx By Bz Tx Ty Tz mx my mz Ex Ey Ez loss eJx eJy eJz  pBz"<< endl;}
	f << "#step: " << it-1 << endl;

	for (k=0;k<ncz;k++){
			for (j=0;j<ncy;j++){
					for (i=0;i<ncx;i++){
								m=icell(i, j, k);	

							  vec_cpy(J, vc[m].J);
							  Jc1 = vc[m].Jc;
								Jcpar = vc[m].Jcpar;
								Nfac = vc[m].N;
							  vec_cpy(A, vc[m].A);
							  vec_cpy(B_J, vc[m].B_J);			
							  vec_cpy(B, vc[m].B);		
								vec_cpy(T, vc[m].T);
								vec_cpy(mm, vc[m].m);
								vec_cpy(E, vc[m].E);

							  U=vc[m].U;
								P=vc[m].loss;
								vec_cpy(e,vc[m].uv);

								f<< m <<"\t"<< vc[m].rc[0] <<"\t"<< vc[m].rc[1] <<"\t"<< vc[m].rc[2] <<"\t"<< vc[m].rc[3] <<"\t"<< vc[m].rc[4] <<"\t"<< vc[m].rc[5] <<"\t"<< Jcpar 
											<<"\t"<< J[0] <<"\t"<< J[1] <<"\t"<< J[2] <<"\t"<< vec_mod(J) <<"\t"<< Jc1 << "\t" << U <<"\t"<< Nfac << "\t"<< A[0] << "\t" << A[1]<<"\t"<< A[2] 
											<<"\t"<< B_J[0] << "\t" << B_J[1]<<"\t"<< B_J[2] <<"\t"<< B[0] << "\t" << B[1]<<"\t"<< B[2] <<"\t"<< T[0] << "\t" << T[1]<<"\t"<< T[2] 
											<<"\t"<< mm[0] << "\t" << mm[1]<<"\t"<< mm[2] <<"\t"<< E[0] <<"\t"<< E[1] <<"\t"<< E[2] <<"\t"<< P <<"\t"<< e[0] << "\t" << e[1]<<"\t"<< e[2] <<"\t"<< vc[m].pB[2] << endl;				
/*
								f<< m <<"\t"<< vc[m].rc[0] <<"\t"<< vc[m].rc[1] <<"\t"<< vc[m].rc[2] <<"\t"<< vc[m].rc[3] <<"\t"<< vc[m].rc[4] <<"\t"<< vc[m].rc[5] <<"\t"<< Jcpar 
											<<"\t"<< J[0] <<"\t"<< J[1] <<"\t"<< J[2] <<"\t" << endl;				
*/
							}}
											f<< endl;
											f<< endl;}
	f.close();

	f.open("data/cross_sectionX.txt",ios::out|ios::app);

	n=0;

	if (it==1) {f <<"#1 2     3     4 		5  6  7  8  9"<< endl;
							f <<"#n rc[x] rc[y] rc[z] Jx Ex Ey Ez P"<< endl;}
	f << "#step: " << it-1 << endl;

	for (k=0;k<ncz;k++){
			for (j=0;j<ncy;j++){
					for (i=0;i<ncx;i++){
								m=icell(i, j, k);	
							  vec_cpy(J, vc[m].J);
								vec_cpy(E, vc[m].E);
								P=vc[m].loss;

								if (i==(ncx/2)){f<< n <<"\t"<< vc[m].rc[0] <<"\t"<< vc[m].rc[1] <<"\t"<< vc[m].rc[2] <<"\t"<< J[0] <<"\t"<< E[0] <<"\t"<< E[1] <<"\t"<< E[2]<<"\t"<<P<<endl;				
																	n++;}}}
											f<< endl;
											f<< endl;}
	f.close();

	n=0;I=0;

	f.open("data/cross_sectionY.txt",ios::out|ios::app);

	if (it==1) {f <<"#1 2     3     4 		5  6  7  8  9 10 "<< endl;
							f <<"#n rc[x] rc[y] rc[z] Jy Ex Ey Ez P Bz"<< endl;}
	f << "#step: " << it-1 << endl;

	for (k=0;k<ncz;k++){
			for (j=0;j<ncy;j++){
					for (i=0;i<ncx;i++){
								m=icell(i, j, k);	
							  vec_cpy(J, vc[m].J);
								vec_cpy(E, vc[m].E);
							  vec_cpy(B, vc[m].B);
								P=vc[m].loss;

								if (j==(ncy/2)){f<< n <<"\t"<< vc[m].rc[0] <<"\t"<< vc[m].rc[1] <<"\t"<< vc[m].rc[2] <<"\t"<< J[1] <<"\t"<< E[0] <<"\t"<< E[1] <<"\t"<< E[2]<<"\t"<<P<<"\t"<<B[2]<<endl;				
																n++;	I+=vc[m].J[1] * vc[m].size[0] * vc[m].size[2];}}}}
	f<< endl;
	f<< endl; 
	f.close();

	cout << "Total current " << I << " A " << endl;

	f.open("data/output3Dx.txt",ios::out|ios::app);

	if(it==1){f <<"#1 2     3     4  		5  		 6   		 7   		8   	 9  10 11 12  13  14  15 	16   17   18 	 19 	20 	 21 	22 23 24 25 26 27 28 29 30 31 32 33	34	 35	 36	 37"<< endl;
						f <<"#n rc[x] rc[y] rc[z] rec[r] rec[fi] rec[z] Jcpar  Jx Jy Jz |J| Jc  U   N 	A_ax A_ay A_az B_Jx B_Jy B_Jz Bx By Bz Tx Ty Tz mx my mz Ex Ey Ez loss eJx eJy eJz"<< endl;}

	f << "#step: " << it-1 << endl;

	for (i=0;i<ncx;i++){
			for (k=0;k<ncz;k++){
					for (j=0;j<ncy;j++){
								m=icell(i, j, k);	

							  vec_cpy(J, vc[m].J);
							  Jc1 = vc[m].Jc;
							  vec_cpy(A, vc[m].A);
							  vec_cpy(B_J, vc[m].B_J);
							  vec_cpy(B, vc[m].B);
								vec_cpy(T, vc[m].T);
								vec_cpy(mm, vc[m].m);				
								P=vc[m].loss;
								vec_cpy(e,vc[m].uv);					

								f<< m <<"\t"<< vc[m].rc[0] <<"\t"<< vc[m].rc[1] <<"\t"<< vc[m].rc[2] <<"\t"<< vc[m].rc[3] <<"\t"<< vc[m].rc[4] <<"\t"<< vc[m].rc[5] <<"\t"<< Jcpar
											<<"\t"<< J[0] <<"\t"<< J[1] <<"\t"<< J[2] <<"\t"<< vec_mod(J) <<"\t"<< Jc1 <<"\t"<< vc[m].U <<"\t"<< vc[m].N <<"\t"<< A[0] <<"\t" << A[1] <<"\t"<< A[2] //<< endl; 
											<<"\t"<< B_J[0] <<"\t"<< B_J[1] <<"\t"<< B_J[2] <<"\t"<< B[0] <<"\t"<< B[1] <<"\t"<< B[2] <<"\t"<< T[0] <<"\t"<< T[1] <<"\t"<< T[2] <<endl;
	//										<<"\t"<< mm[0] <<"\t"<< mm[1] <<"\t"<< mm[2] <<"\t"<< E[0] <<"\t"<< E[1] <<"\t"<< E[2] <<"\t"<< P <<"\t"<< e[0] <<"\t" << e[1] <<"\t"<< e[2] <<endl;						
	/*
								f<< m <<"\t"<< vc[m].rc[0] <<"\t"<< vc[m].rc[1] <<"\t"<< vc[m].rc[2] <<"\t"<< vc[m].rc[3] <<"\t"<< vc[m].rc[4] <<"\t"<< vc[m].rc[5] <<"\t"<< Jcpar
											<<"\t"<< J[0] <<"\t"<< J[1] <<"\t"<< J[2] <<"\t"<<endl;						
	*/
													}}
									f<< endl;
									f<< endl; 
									}
	f.close();

	f.open("data/output3Dy.txt",ios::out|ios::app);

	if (it==1) {f <<"#1 2     3     4  		5  	 	 6   		 7   		8   	9  10 11 12  13  14  15	 16   17   18 	19 	 20 	21 	 22 23 24 25 26 27 28 29 30 31 32 33 34	  35	36	37"<< endl;
							f <<"#n rc[x] rc[z] rc[y] rec[r] rec[fi] rec[z] Jcpar Jx Jy Jz |J| Jc  0	 0 	 A_ax A_ay A_az B_Jx B_Jy B_Jz Bx By Bz Tx Ty Tz mx my mz Ex Ey Ez loss eJx eJy eJz"<< endl;}
	f << "#step: " << it-1 << endl;

	for (j=0;j<ncy;j++){
			for (k=0;k<ncz;k++){
					for (i=0;i<ncx;i++){
								m=icell(i, j, k);	

							  vec_cpy(J, vc[m].J);
							  Jc1= vc[m].Jc;
							  vec_cpy(A, vc[m].A);	
							  vec_cpy(B_J, vc[m].B_J);	
							  vec_cpy(B, vc[m].B);	
								vec_cpy(T, vc[m].T);
								vec_cpy(mm, vc[m].m);	
								P=vc[m].loss;
								vec_cpy(e,vc[m].uv);
							
								f<< m <<"\t"<< vc[m].rc[0] <<"\t"<< vc[m].rc[1] <<"\t"<< vc[m].rc[2] <<"\t"<< vc[m].rc[3] <<"\t"<< vc[m].rc[4] <<"\t"<< vc[m].rc[5] <<"\t"<< Jcpar
											<<"\t"<< J[0] <<"\t"<< J[1] <<"\t"<< J[2] <<"\t"<< vec_mod(J) <<"\t"<< Jc1 <<"\t"<< 0	<<"\t"<< 0 <<"\t"<< A[0] <<"\t"<< A[1]<<"\t"<< A[2] 
											<<"\t"<< B_J[0] <<"\t"<< B_J[1] <<"\t"<< B_J[2] <<"\t"<< B[0] <<"\t"<< B[1]<<"\t"<< B[2]<<"\t"<< T[0] <<"\t"<< T[1]<<"\t"<< T[2] 
											<<"\t"<< mm[0] <<"\t"<< mm[1] <<"\t"<< mm[2] <<"\t"<< E[0] <<"\t"<< E[1] <<"\t"<< E[2] <<"\t"<< P <<"\t"<< e[0] <<"\t"<< e[1] <<"\t"<< e[2] <<endl;
	/*
								f<< m <<"\t"<< vc[m].rc[0] <<"\t"<< vc[m].rc[1] <<"\t"<< vc[m].rc[2] <<"\t"<< vc[m].rc[3] <<"\t"<< vc[m].rc[4] <<"\t"<< vc[m].rc[5] <<"\t"<< Jcpar
											<<"\t"<< J[0] <<"\t"<< J[1] <<"\t"<< J[2] <<endl;
		*/						
															}}
											f<< endl;
											f<< endl; 
											}
	f.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_plane()
{int i,j,k,m, n=0;
 double J[3]={0}, Jc[3]={0},Na=0, A[3]={0}, B[3]={0}, B_J[3]={0}, T[3]={0}, mm[3]={0}, P=0, e[3]={0},E[3]={0},U=0;

	fstream f;
	f.open("data/output_plane.txt",ios::out|ios::app);
	if (it==1) {f <<"#1 2     3     4 		5   6   7   8  9  10"<< endl;
							f <<"#n rc[x] rc[y] rc[z] BJx BJy BJz Bx By Bz"<< endl;}
	f << "#step: " << it-1 << endl;

	for (k=0;k<ncz_plane;k++){
			for (j=0;j<ncy_plane;j++){
					for (i=0;i<ncx_plane;i++){

								m=iplanecell(i, j, k);	
								f<< m <<"\t"<< RC_plane[m][0] <<"\t"<< RC_plane[m][1] <<"\t"<< RC_plane[m][2] <<"\t"<< BJ_plane[m][0] <<"\t"<< BJ_plane[m][1] <<"\t"<< BJ_plane[m][2] << 
												"\t"<< B_plane[m][0] <<"\t"<< B_plane[m][1] <<"\t"<< B_plane[m][2] <<endl;}}
								f<< endl;
								f<< endl;}
	f.close();

	f.open("data/cross_plane.txt",ios::out|ios::app);

	n=0;

	if (it==1) {f <<"#1 2 3     4 		5     6"<< endl;
							f <<"#n t rc[x] rc[y] rc[z] Bz"<< endl;}
	f << "#step: " << it-1 << endl;

	for (k=0;k<ncz_plane;k++){
			for (j=0;j<ncy_plane;j++){
					for (i=0;i<ncx_plane;i++){

								m=iplanecell(i, j, k);	
								if(j==ncy/2){
											f<< m <<"\t"<< t <<"\t"<< RC_plane[m][0] <<"\t"<< RC_plane[m][1] <<"\t"<< RC_plane[m][2] <<"\t"<< BJ_plane[m][0] <<"\t"<< BJ_plane[m][1] <<"\t"<< BJ_plane[m][2] << 
															"\t"<< B_plane[m][0] <<"\t"<< B_plane[m][1] <<"\t"<< B_plane[m][2] <<endl;}}}
								f<< endl;
								f<< endl;}
	f.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_surface()
{int n,i,j,k;

	fstream f;
	f.open("data/output_Jx.txt",ios::out|ios::app);
	if (it==1) f << "#1	2	3	4	5		 6		7		 8 9	10"<< endl;
	if (it==1) f << "#n i j k X[0] X[1] X[2] J dJ av_dA"<< endl;
	f << "#step: " << it-1 << endl;

	for (k=0;k<ncz;k++){
			for (j=0;j<ncy;j++){
					for (i=0;i<=ncx;i++){
								n=ijktonsurX(i,j,k);
								f<<n<<"\t"<< i <<"\t"<< j <<"\t"<< k <<"\t"<<vsX[n].rc[0]<<"\t"<<vsX[n].rc[1]<<"\t"<<vsX[n].rc[2]<<"\t"<<vsX[n].Jo + vsX[n].dJ <<"\t"<< vsX[n].dJ<<"\t"<<vsX[n].av_dA<<endl;}}
//			f<< endl;
//			f<< endl;
}
	f.close();

	f.open("data/output_Jy.txt",ios::out|ios::app);
	if (it==1) f << "#1	2	3 4 5    6    7    8 9	10		 11  12   13" << endl;	
	if (it==1) f << "#n i j k X[0] X[1] X[2] J dJ av_dAa Aa0 snAy us"<< endl;
	f << "#step: " << it-1 << endl;

	for (k=0;k<ncz;k++){
			for (j=0;j<=ncy;j++){
					for (i=0;i<ncx;i++){
								n=ijktonsurY(i,j,k);
								if(Ismax!=0){f<<n<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<vsY[n].rc[0]<<"\t"<<vsY[n].rc[1]<<"\t"<<vsY[n].rc[2]<<"\t"<<vsY[n].Jo+vsY[n].dJ<<"\t"<<vsY[n].dJ<<"\t"<<vsY[n].av_dA<<"\t"<<vsY[n].Aa0 <<"\t"<<snAy[n]<<"\t"<<vsY[n].us<< endl;}
								else{f<<n<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<vsY[n].rc[0]<<"\t"<<vsY[n].rc[1]<<"\t"<<vsY[n].rc[2]<<"\t"<<vsY[n].Jo+vsY[n].dJ<<"\t"<<vsY[n].dJ<<"\t"<<vsY[n].dAa<<"\t"<< vsY[n].Aa0 << endl;}}}
//			f<< endl;
//			f<< endl;
}
	f.close();

	f.open("data/output_Jz.txt",ios::out|ios::app);
	if (it==1) f << "#1	2	3 4 5    6    7    8 9  10    11"<< endl;  	
	if (it==1) f << "#n i j k X[0] X[1] X[2] J dJ av_dA snAz"<< endl;
	f << "#step: " << it-1 << endl;
	for (k=0;k<=ncz;k++){
			for (j=0;j<ncy;j++){
					for (i=0;i<ncx;i++){
								n=ijktonsurZ(i,j,k);
								if(Ismax!=0){f<<n<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<vsZ[n].rc[0]<<"\t"<<vsZ[n].rc[1]<<"\t"<<vsZ[n].rc[2]<<"\t"<<vsZ[n].Jo + vsZ[n].dJ <<"\t"<<vsZ[n].dJ<<"\t"<<vsZ[n].av_dA<< "\t"<<snAz[n]<<endl;}
								else{f<<n<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<vsZ[n].rc[0]<<"\t"<<vsZ[n].rc[1]<<"\t"<<vsZ[n].rc[2]<<"\t"<<vsZ[n].Jo + vsZ[n].dJ <<"\t"<<vsZ[n].dJ<<"\t"<<vsZ[n].av_dA<< endl;}
				}}
//			f<< endl;
//			f<< endl;
				}
	f.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::cal_AV()
{int i,j,k,n,m;

		for(j=0;j<ncy;j++){
				for(i=0;i<ncx;i++){
										n=icell(i, j, 0);
										Jx[n]=0;Jy[n]=0;Jz[n]=0;
										Bx[n]=0;By[n]=0;Bz[n]=0;
										Tx[n]=0;Ty[n]=0;Tz[n]=0;}}


		for(k=0;k<ncz;k++){
				for(j=0;j<ncy;j++){
						for(i=0;i<ncx;i++){
										m=icell(i, j, k);	
										n=icell(i, j, 0);

										Jx[n]+=vc[m].J[0]/ncz;			
										Jy[n]+=vc[m].J[1]/ncz;
										Jz[n]+=vc[m].J[2]/ncz;
										Bx[n]+=vc[m].B[0]/ncz;
										By[n]+=vc[m].B[1]/ncz;
										Bz[n]+=vc[m].B[2]/ncz;
										Tx[n]+=vc[m].T[0]/ncz;
										Ty[n]+=vc[m].T[1]/ncz;
										Tz[n]+=vc[m].T[2]/ncz;}}}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_AV()	
{int i,j,k,n,m;

	fstream f;
	f.open("data/outputAV.txt",ios::out|ios::app);
	if (it==1) {f <<"#1 2     3     4 		5  6  7  8  9  10 11 12 13"<< endl;
							f <<"#n rc[x] rc[y] rc[z] Jx Jy Jz Bx By Bz Tx Ty Tz"<< endl;}
	f << "#step: " << it-1 << endl;

			k=0;
			for(j=0;j<ncy;j++){
					for(i=0;i<ncx;i++){
								m=icell(i, j, k);	

								f<< m <<"\t"<< vc[m].rc[0] <<"\t"<< vc[m].rc[1] <<"\t"<< vc[m].rc[2] <<"\t"<< Jx[m] <<"\t"<< Jy[m] <<"\t"<< Jz[m] <<"\t"<< Bx[m] <<"\t"<< By[m] <<"\t"<< Bz[m] <<
												"\t"<< Tx[m] <<"\t"<< Ty[m] <<"\t"<< Tz[m] <<endl;}}	
						f<< endl;
						f<< endl; 
						f.close();

	k=0;
	f.open("data/Tz.txt",ios::out|ios::app);{	
			f << "#step: " << it-1 << endl;
			for(i=0;i<ncx;i++){
						for(j=0;j<ncy;j++){		
										n=icell(i,j,k);
										f<< vc[n].rc[0] << "\t" << vc[n].rc[1] << "\t" << Tz[n] << endl;}	
						f << endl;}

			f << endl;

			for(j=0;j<ncy;j++){
						for(i=0;i<ncx;i++){		
										n=icell(i,j,k);
										f<< vc[n].rc[0] << "\t" << vc[n].rc[1] << "\t" << Tz[n] << endl;}	
						f << endl;}
						f << endl;
						f << endl;}
	f.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void class_cube::save_Tz()	
{int i,j,k,n,m;

	k=0;
	fstream f;
	f.open("data/Tz.txt",ios::out|ios::app);{	
			f << "#step: " << it-1 << endl;
			for(i=0;i<ncx;i++){
						for(j=0;j<ncy;j++){		
										n=icell(i,j,k);
										f<< vc[n].rc[0] << "\t" << vc[n].rc[1] << "\t" << vc[n].T[2] << endl;}	
								f << endl;}

			f << endl;

			for(j=0;j<ncy;j++){
						for(i=0;i<ncx;i++){		
										n=icell(i,j,k);
										f<< vc[n].rc[0] << "\t" << vc[n].rc[1] << "\t" << vc[n].T[2] << endl;}	
										f << endl;}
						f << endl;
						f << endl;}

			f.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
