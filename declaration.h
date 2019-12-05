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

struct struct_node			
{	
	int adr[3];												//i,j,k address 
	double rc[6];											//vector of nodes x,y,z,r,fi,z
}; 	

struct struct_edgeX
{
	int adr[3];												//i,j,k address 
	int adrs[4];											//adresses of surfaces y1,y2,z1,z2
	int adrc[4];											//adresses of cells m1,m2,m3,m4
	int adres[7];											//adresses of symmetry edges X 
	double rc[6];											//vector position of edge RC[m]={x, y, z, r, FI, z}
	double l;													//lenght of the edge [m]
	double T0;												//Tx component
	double dT;												//Tx component
	double dTx;												//temporary variable
};

struct struct_edgeY
{
	int adr[3];												//i,j,k address 
	int adrs[4];											//adresses of surfaces x1,x2,z1,z2
	int adrc[4];											//adresses of cells m1,m2,m3,m4
	int adres[7];											//adresses of symmetry edges Y
	double rc[6];											//vector position of edge RC[m]={x, y, z, r, FI, z}
	double l;													//lenght of the edge [m]
	double T0;												//Ty component
	double dT;												//Ty component
	double dTx;												//temporary variable
};

struct struct_edgeZ
{
	int adr[3];												//i,j,k address 
	int adrs[4];											//adresses of surfaces x1,x2,y1,y2
	int adrc[4];											//adresses of cells m1,m2,m3,m4
	int adres[7];											//adresses of symmetry edges Z 
	double rc[6];											//vector position of edge RC[m]={x, y, z, r, FI, z}
	double l;													//lenght of the edge [m]
	double T0;												//Tz component
	double dT;												//Tz component
	double dTx;												//temporary variable
};

struct struct_surfaceX
{
	int adr[3];												//i,j,k address 
	int adre[4];											//adresses of edges y1,y2,z1,z2
	int adrss[7];											//adresses of symmetry surfaces X
	int sadr[3];											//address of sector for multipole expansion
	int rijk[6];											//range of surfaces for inner boxes of multipole expansion 
	int na;														//total number of surfaces for inner boxes of the multipole expansion
	int *n;														//surface addreesses for inner boxes of the multipole expansion
	double *r3;												//surface distance to the sector for the multipole expansion
	double *rx;												//surface x distance to the sector for the multipole expansion
	double *ry;												//surface y distance to the sector for the multipole expansion
	double *rz;												//surface z distance to the sector for the multipole expansion
	double rc[6];						  				//vector position of surface RC[m]={x, y, z, r, FI, z}
	double size[4];										//size of cell SS[m]={a, b, c} in polar SS[m]={dr, rdfi, dz, dfi}
	double S;													//area of surface X
	double Jo;												//current density at time t_o Jo[A]={jx, jy, jz}  
	double dJ;		   									//delta current density dJ[A]={jx, jy, jz} total current density: J = Jo + dJ
	double dJx;												//temporary variable for current density
	double J;													//total current density
	double Vi;	      								//volume of influence
	double Aa0;	  										//magnetic vector potential Ao[T.m]=ax 
	double dAa;												//delta magnetic vector potential dA[T.m]=ax 
	double av_dA;											//delta average potential
	double dAx;							 					//temporary variable delta potential 
	double us;												//unit vector of source current
	double Is;												//source current
	double dI;												//delta of source current
	double dAs;												//delta vector potential of source current
	int tmp;
};

struct struct_surfaceY
{
	int adr[3];												//i,j,k address 
	int adre[4];											//adresses of edges x1,x2,y1,y2
	int adrss[7];											//adresses of symmetry surfaces Y
	int sadr[3];											//address of sector for multipole expansion
	int rijk[6];											//range of surfaces for inner boxes of multipole expansion 
	int na;														//total number of surfaces for inner boxes of the multipole expansion
	int *n;														//surface addreesses for inner boxes of the multipole expansion
	double *r3;												//surface distance to the sector for the multipole expansion
	double *rx;												//surface x distance to the sector for the multipole expansion
	double *ry;												//surface y distance to the sector for the multipole expansion
	double *rz;												//surface z distance to the sector for the multipole expansion
	double rc[6];								  		//vector position of surface RC[m]={x, y, z, r, FI, z}
	double size[4];										//size of cell SS[m]={a, b, c} in polar SS[m]={dr, rdfi, dz, dfi} 
	double S;													//area of surface Y
	double Jo;												//current density at time t_o Jo[A]={jx, jy, jz}  
	double dJ;		   									//delta current density dJ[A]={jx, jy, jz} total current density: J = Jo + dJ
	double dJx;												//temporary variable for current density
	double J;													//total current density
	double Vi;												//volume of influence
	double Aa0;	  										//magnetic vector potential Ao[T.m]=ay
	double dAa;												//delta magnetic vector potential dA[T.m]=ay 
	double av_dA;											//delta average potential
	double dAx;							 					//temporary variable delta potential 
	double us;												//unit vector of source current
	double Is;												//source current
	double dI;												//delta of source current
	double dAs;												//delta vector potential of source current
};

struct struct_surfaceZ
{
	int adr[3];												//i,j,k address 
	int adre[4];											//adresses of edges x1,x2,y1,y2
	int adrss[7];											//adresses of symmetry surfaces Z
	int sadr[3];											//address of sector for multipole expansion
	int rijk[6];											//range of surfaces for inner boxes of multipole expansion
	int na;														//total number of surfaces for inner boxes of the multipole expansion
	int *n;														//surface addreesses for inner boxes of the multipole expansion
	double *r3;												//surface distance to the sector for the multipole expansion
	double *rx;												//surface x distance to the sector for the multipole expansion
	double *ry;												//surface y distance to the sector for the multipole expansion
	double *rz;												//surface z distance to the sector for the multipole expansion
	double rc[6];								  		//vector position of surface RC[m]={x, y, z, r, FI, z}
	double size[4];										//size of cell SS[m]={a, b, c} in polar SS[m]={dr, rdfi, dz, dfi}
	double S;													//area of surface Z
	double Jo;												//current density at time to Jo[A/m2]={jx, jy, jz}  
	double dJ;		   									//delta current density dJ[A/m2]={jx, jy, jz} total current density: J = Jo + dJ
	double dJx;												//temporary variable for current density
	double J;													//total current density
	double Vi;												//volume of influence
	double Aa0;	  										//magnetic vector potential Ao[T.m]=az
	double dAa;												//delta magnetic vector potential dA[T.m]=az 
	double av_dA;											//delta average potential
	double dAx;							 					//temporary variable delta potential
	double us;												//unit vector of source current
	double Is;												//source current
	double dI;												//delta of source current
	double dAs;												//delta vector potential of source current
};

struct struct_cell  
{
	int adr[3];												//i,j,k address 
	int adrp[3];											//i,j,k address in polar coordinate system 
	double rc[6];			    						//vector position of cells RC[m]={x, y, z, r, FI, z}
	double size[3];										//size of cell SC[m]={a, b, c}
	double V;				     							//volume of cell V[m]
	double J[3];											//current density J[A]=[jx, jy, jz] after interpolation of surfaces {Jx, Jy, Jz} 
	double uv[3];				              //unit vector of current per cell
	double loss;											//local loss per cell P[W]
	double B_J[3];										//magnetci field B[T]={Bx, By, Bz}	created by current density
	double B[3];											//total magnetic field B (magnetic field + applied field)	
	double theta;											//angel of the magnetic field with z axis
	double pB[3];											//total magnetic field B (magnetic field + applied field)	
	double A[3];											//magnetic vector potential field B[T.m] 	
	double Jc;												//critical current density Jc[A/m2]={Jcx, Jcy, Jcz}
	double Jcpar;											//critical current density Jc[A/m2]={Jcx, Jcy, Jcz}
	double Jcx;												//temporary variable for critical current density
	double N;													//N factor for power law
	double T[3];											//T-vector
	double m[3];											//magnetic moment
	double U;													//U dissipation factor
	double E[3];											//electric field	
	double Js[3];											//current density from current source	
};

class class_cube      
{
	protected:

  //Input variables//	
	double x;													//size of x [m]
	double xl;												//thickness of linear material [m]
	double y;													//size of y [m]
 	double z;													//size of z [m]
	int full_matrix;									//full size interaction matric											
	int nsucx;												//number of superconducting cells x 
	int nncx;													//number of normal (linear) cells x
	int ncx;													//number of cells x 
	int ncy;   		 										//number of cells y 
	int ncz;													//number of cells z 
	int n_tapes;											//number of tapes in stack 
	int nc_tape;											//number of cells in 1 tape 
	int nc_gap;												//number of cells in 1 gap between stack 
	double d_tape;										//thickness of 1 tape in stack 
	double d_gap;											//thickness of gap between two tapes 
	int elc;													//calculation of vetor potential by approximation or numerical way 
	double tol_elc;										//tolerance for auto numerical calculation of average vector potantial 
	double Bamax;											//aplied magnetic field Bamax[T]
	double Bamax1;										//aplied magnetic field Bamax[T] for cross-field emagnetzation in different direction as magnetizate applied field
	int Bshape;												//shape of applied magnetic field
	int Btrape;												//trapped field in the plane outside of sample
	double Ismax;											//current from current source
	double rcx_plane;									//position vector of the plane	
	double rcy_plane;									//position vector of the plane
	double rcz_plane;									//position vector of the plane
	double x_plane;										//size of the plane
	double y_plane;										//size of the plane
	double z_plane;										//size of the plane
	int ncx_plane;										//number of cells along X axis in the plane
	int ncy_plane;										//number of cells along X axis in the plane
	int ncz_plane;										//number of cells along X axis in the plane
	double fi;												//angle in spherical coordinates
	double fi1;												//angle in spherical coordinates for cross-field demagnetization in different direction as magnetizate applied field
	double theta;											//angle in spherical coordinates	
	double uni;												//mesh: 1 - uniform mesh, 0 - non-uniform mesh, 2 - semi-uniform mesh
	double rel;												//relation: 1 - Jo constant; 2 - Jc(B) Kim; 3 - Jc(B) measured data; 4 - anisotropic; 5 - CSM
	int nB;														//0-no,1-yes with n(B) dependence
	int measured_points;							//number of measured points in Jc(B) curve
	int measured_fields;							//number of measured B fields in Jc(B,theta) curve
	double sym;												//symmetric calculation
	double Ec;												//critical electric field Ec[V/m]
	double Jo;												//critical current density Jc[A/m2]
	double Jol;												//critical current density Jc[A/m] for linear material
	double rhoR;											//real resistivity for linear material
	double dl;												//real thickness of linear material	
 	double Jcpa;											//parallel critical current density
	double Jcpe;											//perpendicular critical current density
	double Bo;												//Bo Kim model 20mT		
	double N;													//N order of power law
	double Nl;												//N order of Power law for linear material N=1 (Ohm's law)		
	double m;													//m Kim model 0.5
	double f;													//frequency
	double f1;												//frequency
	int ns;														//number of samples per one period						
	double tola;											//tolerance of small number
	int step;													//number of steps of dt in cycle
	int smooth;												//interpolated current density higher resolution
	int shape;												//shape of thin film
	int num_threads;									//number of threads	
	int shiftx1;											//shift of sectors set2 along Z axes 			
	int shiftx2;											//shift of sectors set3 along Z axes
	int shifty1;											//shift of sectors set2 along Z axes 			
	int shifty2;											//shift of sectors set3 along Z axes
	int shiftz1;											//shift of sectors set2 along Z axes 			
	int shiftz2;											//shift of sectors set3 along Z axes


	int nsx[3];												//numbers of sectors in X direction for three sets 
	int nsy[3];												//numbers of sectors in Y direction for three sets
	int nsz[3];												//numbers of sectors in Z direction for three sets
	int NCX;													//total number of cells in sectors along the X axis
	int NCY;													//total number of cells in sectors along the Y axis
	int NCZ;													//total number of cells in sectors along the Z axis	

	//sectors//
	int set;													//total number of sets set=0 not shifted, set=1 shift1, set=2 shift2	
	int nsall;												//total number of sectors (set0+set1+set2)
	int xsall[3];											//total number of not shifterd sectorsxsall[0]=set0,....
	int nscx;													//total number of X cells in one sector
	int nscy;													//total number of Y cells in one sector
	int nscz;													//total number of Z cells in one sector
  
	int *cii, *cai;										//minimum and maximum of i at cells
	int *cij, *caj;										//minimum and maximum of j at cells
	int *cik, *cak;										//minimum and maximum of k at cells


	//cells//
	double lx; 												//size of cell x[m]
	double ly; 												//size of cell y[m]
	double lz; 												//size of cell z[m]
	double cx; 												//size of cell x[m]
	double cy; 												//size of cell y[m]
	double cz; 												//size of cell z[m]
	int nc;														//total number of all cells
 	int nall;        									//total number of elements of matrix n*n
	int cube;													//shape of sample	

	//surfaces//

	int nsurx;												//total number of surfaces X direction
	int nsury;												//total number of surfaces Y direction	
	int nsurz;												//total number of surfaces Z direction
	int nedgex;												//total number of edge X
	int nedgey;												//total number of edge Y
	int nedgez;												//total number of edge Z
	int nsxall;												//total number of elements of matrix nsurx*nsurx
	int nsyall;												//total number of elements of matrix nsury*nsury
	int nsxyall;											//total number of elements of matrix nsurx*nsury
	int nszall;												//total number of elements of matrix nsurz*nsurz
	int nxall;												//total number of elements of matrix nc*nsurx
	int nyall;												//total number of elements of matrix nc*nsury
	int nzall;												//total number of elements of matrix nc*nsurz

  //nodes//
	int nnx;													//number of x nodes 
	int nny;													//number of y nodes
	int nnz;													//number of z nodes
	int nn;														//total number of all nodes
	
	//constants//
	double q;													//charge
	double pi;												//3.1415	
	double ep0;												//epsilon
	double mi0;												//mi0
	double c;													//speed of light
	
	//time//
	double t;													//time [s]
	double tp;												//time [s] in previous itteration step
	double dt;												//delta time
	double T;													//period of the signal
	double dt1;												//delta time
	double T1;												//period of the signal
	int it;														//number of itteration step
	double Ba;												//instant applied magnetic field [T] 
	double B_a[3];										//instant applied magnetic field [T] vector include angle
	double pBa[3]; 										//previous applied magnetic field [T] 
	double dBa[3];										//delta applied magnetic field [T] 
	double ea[3];											//unit vector of applied field
	double Is;												//current from current source	
	double pI;												//previous time step current from current source	
	double dI;												//delta current from current source	

	//AC loss//
	double Px;												//loss in one dt
	double Pix;												//loss in previous dt
	double Plin;											//loss in linear material in one dt
	double Psup;											//loss in superconducting material in one dt
	double Pilin;											//linear loss in previous dt
	double Pisup;											//superconducitng loss in previous dt
	double Ph;												//loos in one dt thought hole cycle	
	double Q;													//total loss in sampel per entire cycle	
	double Qs;												//total loss in superconducting in sampel per entire cycle	
	double Ql;												//total loss in linear material in sampel per entire cycle	
	double Qh;												//total loss calculated from hysterezis loop 	
	double Qh1;												//total loss calculated from magnetizatioin first half loop 
	double mp[3];											//previous magnetic moment

	//Other//
	double zero[3];										//zero vector
	double ratio; 										//ratio of sample size x:y		(1/(x/y))
	double d_ncx;											//number of cells x 
	double d_ncy;											//number of cells y		
	double l_debug;										//debuging variable				
	double **RC_plane;
	double **BJ_plane;
	double **B_plane;

	//minimization//
	int stepall;											//number of steps of minimization loop
	double tolJ;											//tolerance of currnet density 
	double hp;								  			//constant of pozitive dJx 
	double hn;								  			//constant of negative dJx
	double htpz;											//constant of pozitive dT		
	double htnz;  										//constant of negative dT
	double htpy;											//constant of pozitive dT		
	double htny;  										//constant of negative dT
	double htpx;											//constant of pozitive dT		
	double htnx;  										//constant of negative dT
	double V_ave;											//average volume
	double V_tot;											//total volume
	double minJc;											//minimum value of critical current density
	double Uo;												//U constant for anisotropic case
	double mo;                        //m constant for anisotropic case
	double rho;											  //resistivity	

	//flux line//									
	double hr;												//distance between two steps of one flux line

	//plane//
	int	nc_plane; 										//total number of cells on the plane

	double plx;												//size of cell x[m]
	double ply;												//size of cell y[m]
	double plz;												//size of cell z[m]
	double pcx[3];										//size of cell x[m] in plane pcx[x,y,z]
	double V_ave_plane;								//average volume of cell in plane

	//matrix in memory//
	struct_node *vnode;								//variabels of nodes
	struct_edgeX *veX;								//variables of edges X
	struct_edgeY *veY;								//variables of edges Y
	struct_edgeZ *veZ;								//variables of edges Z
	struct_surfaceX *vsX;							//variabels of surfaces X
	struct_surfaceY *vsY;							//variabels of surfaces Y
	struct_surfaceZ *vsZ;							//variabels of surfaces	Z
	struct_cell *vc;									//variabels of cells

	double *cmAx;											//potencial matrix for surfaces X size 1 x nsurx
	double *cmAy;											//potencial matrix for surfaces Y size 1 x nsury
	double *cmAz;											//potencial matrix for surfaces Z size 1 x nsurz

	double *snAx;											//vector potencial matrix for surfaces X size 1 x nsurx
	double *snAy;											//vector potencial matrix for surfaces Y size 1 x nsury
	double *snAz;											//vector potencial matrix for surfaces Z size 1 x nsurz

	double **mAx;											//potencial matrix for surfaces X	r
	double **mAy;											//potencial matrix for surfaces Y fi
	double **mAz;											//potencial matrix for surfaces Z z

	double *cmHxyz;										//Hx field matrix for surfaces X size 1 x nsurx
	double *cmHxzy;										//Hx field matrix for surfaces X size 1 x nsurx
	double *cmHyzx;										//Hy field matrix for surfaces Y size 1 x nsury
	double *cmHyxz;										//Hy field matrix for surfaces Y size 1 x nsury
	double *cmHzxy;										//Hz field matrix for surfaces Z size 1 x nsurz
	double *cmHzyx;										//Hz field matrix for surfaces Z size 1 x nsurz

	double **mHxyz;										//Hx field matrix for surfaces X size 1 x nsurx
	double **mHxzy;										//Hx field matrix for surfaces X size 1 x nsurx
	double **mHyzx;										//Hy field matrix for surfaces Y size 1 x nsury
	double **mHyxz;										//Hy field matrix for surfaces Y size 1 x nsury
	double **mHzxy;										//Hz field matrix for surfaces Z size 1 x nsurz
	double **mHzyx;										//Hz field matrix for surfaces Z size 1 x nsurz

	double *Jx;												//matrix of average Jx
	double *Jy;												//matrix of average Jy
	double *Jz;												//matrix of average Jz

	double *Bx;												//matrix of average Bx
	double *By;												//matrix of average By
	double *Bz;												//matrix of average Bz

	double *Tx;												//matrix of average Tx
	double *Ty;												//matrix of average Ty
	double *Tz;												//matrix of average Tz

	double **emHxyz;
	double **emHxzy;
	double **emHyzx;
	double **emHyxz;
	double **emHzxy;
	double **emHzyx;

	double **data1;
	double ***data;

	public:
///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Geometry //////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

void input();																						//set input variables 

  ///Nodes///
void node_rc(int i, int j, int k, double XN[3]); 				//vector position of node rc[3] from i,j,k

int  inode(int i,int j,int k);				         					//address of node n from i,j,k
void nnodetoijk(int n, int& i, int& j, int& k);         //address of node i,j,k from n

  ///Edges///
void edgeX_rc(int i, int j, int k, double RC[3]);				//vector position of the edge rc[3]] in the direction x 
void edgeY_rc(int i, int j, int k, double RC[3]);   		//vector position of the edge rc[3]] in the direction y
void edgeZ_rc(int i, int j, int k, double RC[3]);				//vector position of the edge rc[3]] in the direction z 
double edgeX_lenght(int i, int j, int k);								//calculates lenght of edge X [m] 
double edgeY_lenght(int i, int j, int k);								//calculates lenght of edge Y [m] 
double edgeZ_lenght(int i, int j, int k);								//calculates lenght of edge Z [m] 

int ijktonedgeX(int i,int j,int k);											//address of edge X n from i,j,k
int ijktonedgeY(int i,int j,int k);											//address of edge Y n from i,j,k
int ijktonedgeZ(int i,int j,int k);											//address of edge Z n from i,j,k
void nedgeXtoijk(int n, int& ix, int& iy, int& iz);			//address of edge X i,j,k from n
void nedgeYtoijk(int n, int& ix, int& iy, int& iz);			//address of edge Y i,j,k from n
void nedgeZtoijk(int n, int& ix, int& iy, int& iz);			//address of edge Z i,j,k from n
void edgeZtosurXY(int m, int& x1, int& x2, int& y1, int& y2);																									 //find 4 neighborhood surfaces to Tz edge
void edgeXtosurYZ(int m, int& y1, int& y2, int& z1, int& z2); 																								 //find 4 neighborhood surfaces to Tx edge
void edgeYtosurXZ(int m, int& x1, int& x2, int& z1, int& z2);																									 //find 4 neighborhood surfaces to Ty edge
void boundary_addres_edge();														//calculates adreesses of cells and surfaces around the edges at the boundary for ring

  ///Surface///
void surfacex_rc(int i, int j, int k, double XSC[3]);		//vector position of the surface X
void surfacey_rc(int i, int j, int k, double YSC[3]);   //vector position of the surface Y
void surfacez_rc(int i, int j, int k, double ZSC[3]);		//vector position of the surface Z
void surfaceX_size (int i, int j, int k, double SS[2]); //size of surface X {a, b} 
void surfaceY_size (int i, int j, int k, double SS[2]);	//size of surface Y {a, b} 
void surfaceZ_size (int i, int j, int k, double SS[2]);	//size of surface Z {a, b}
double surX_S (int m, int i, int j, int k);							//calculates area of surface X					 
double surY_S (int m, int i, int j ,int k);							//calculates area of surface Y	
double surZ_S (int m, int i, int j, int k);							//calculates area of surface Z	
double vol_infX(int i, int j, int k);										//volume of influence of surfaces X
double vol_infY(int i, int j, int k);										//volume of influence of surfaces Y
double vol_infZ(int i, int j, int k);										//volume of influence of surfaces Z
void surf_shape();																			//change shape of sample (disk, ball)

int ijktonsurX(int i,int j,int k);											//address of surface X n from i,j,k
int ijktonsurY(int i,int j,int k);											//address of surface Y n from i,j,k
int ijktonsurZ(int i,int j,int k);											//address of surface Z n from i,j,k
void nsurXtoijk(int n, int& ix, int& iy, int& iz);			//address of surface X i,j,k from n
void nsurYtoijk(int n, int& ix, int& iy, int& iz);			//address of surface Y i,j,k from n
void nsurZtoijk(int n, int& ix, int& iy, int& iz);			//address of surface Z i,j,k from n	

  ///Cells///
double xnode_rc(int i, int j, int k, int l);   					//vector position of nodes according coordinates ijk
double xnode_prc(int i, int j, int k, int l);   				//vector position of nodes according coordinates ijk
void cell_rc (int i, int j, int k, double XC[3]);		 		//vector position of cell 
double cell_rc_plane(int i, int j, int k,int s);				//vector position of cells in the external plane
double cell_vol (double SC[3], int i, int j, int k);		//volume of cell
double cell_size (int i, int j, int k, double SC[3]);		//size of cell 
void cell_lin();																				//set Jol, N in linear cell
int icell(int i,int j,int k);			    		 							//address of cell n from i,j,k
int iplanecell(int i,int j,int k);			    		 				//address of cell n from i,j,k
void ncelltoijk(int n, int& i, int& j, int& k);         //address of cell i,j,k from n
void ncellplanetoijk(int n,int& i,int& j,int& k); 			//address of cell i,j,k from n

  ///Sectors///
void nrange(int nsec);																	//adresses of edge X,Y,Z, cell and surfaces X,Y,Z in n sector
void nrange_db(int nsec);																//cout adresses of cells in all sectors
void cal_sectors();																			//set optimum cells per sector (nscx...) and calculates number of sectors (nsx....)
void nsectortoijk(int s,int n,int& i,int& j,int& k); 		//address of sector i,j,k from n

	/// Interoplation ///
double interX(double r[3], int n);											//interpolation of any point rc from the surface X n  
double interY(double r[3], int n);											//interpolation of any point rc from the surface Y n  
double interZ(double r[3], int n);											//interpolation of any point rc from the surface Z n 
void interJ(double r[3], double J[3]);									//interpolated dJ of any point rc from surfaces   
void interJ_cell(int n, double J[3]);										//interpolated J of cell n (Jo+dJ) 
void interJ0_cell(int n, double J[3]);									//interpolated J0 of cell n Jo 
double interJs_cell(int n, int s);											//interpolated Js of cell n 
void intercell_A(int n, double A[3]);										//interpolated A of cell n
void intercell_T(int n, double T[3]);										//interpolated T of cell n
void interpolation_cell();															//interpolated J, m, A, T, P at the cells
void polar_cart_X(double rc[3], double Jp[3], double J[3]);//trasform  any vctor X, X(r,fi,z) to X(x,y,z)
double bilinear_inter_JcB(double B, double theta);			//linear interpolation of the Jc from measured Jc(B) data
double bilinear_inter_nB(double B, double theta);		  	//linear interpolation of the n from measured n(B) data

	/// Saving data to memory ///
void save_nodes();																			//save node's variables (RC, AN) to memory and file
void save_edges();																			//save edge's variabes (RC, l, T) to memory and file
void save_surfaces();																		//save surface's variables (Vi, RC, size, S, nullified: Jo,dJ,Ao,dA,av_dA) to memory and file
void save_cells ();						 	 												//save cell's variables (size , V, vertex, XC, vol, SC, J) to memory and file

///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Electrodynamic  ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector potential ///
double surX_dA (int i); 																//dA and A from applied field at surface X
double surY_dA (int i); 																//dA and A from applied field at surface Y
double surZ_dA (int i);							 										//dA and A from applied field at surface Z
void cal_av_dA (int s);		 															//average potential at surfaces
void fast_matrix_cmA();																	//average vector potencial matrix of surface X, Y, Z in cartezian coordinate system
void full_matrix_mA();																	//average vector potencial matrix of surface X, Y, Z in cartezian coordinate system
double read_cmAx(int j, int i);													//read from cmAx (vector potential matrix)		
double read_cmAy(int j, int i);													//read from cmAy (vector potential matrix)
double read_cmAz(int j, int i);													//read from cmAz (vector potential matrix)
double num_a_auto(int l, int m, int i); 								//automaticly find the best resolution of sub-elemets
double num_a(int l, int m, int na[3], int nb[3],int s);	//average scalar potential of volume created by volume of 1 element, numerical calculation
void matrix_snA();																			//average vector potencial for source current
void boundary_cmA();																		//average vector potencial matrix of surface X,Y in cartezian coordinate system

	/// Magnetic field ///
void magnetic_field();																	//recalculate magnetic field B= B_J(J) + Ba
void external_magnetic_field();													//recalculate external magnetic field B= B_J(J) + Ba
void magnetic_moment(int n, double m[3]);								//magnetic moment at the cell
double fast_matrix_Hx();																//average magnetic field matrices cmHx, cmHy, cmHz
double full_matrix_Hx();																//average magnetic field matrices cmHx, cmHy, cmHz
double matrix_eHx();																		//average external magnetic field matrices ecmHx, ecmHy, ecmHz		
double read_cmHxyz(int i, int j);												//read from cmHxyz 
double read_cmHxzy(int i, int j);												//read from cmHxzy
double read_cmHyzx(int i, int j);												//read from cmHyzx
double read_cmHyxz(int i, int j);												//read from cmHyxz
double read_cmHzxy(int i, int j);												//read from cmHzxy
double read_cmHzyx(int i, int j);												//read from cmHzxy
void recal_theta();																			//recalculates unit vector of magnetic field in the cells
double cell_theta(int i);																//recalculates unit vector of magnetic field in the cell i
double num_b(int m, int l, int na[3], int s, int t);		//magnetic field of volume created by volume of 1 element, numerical calculation


 /// Current density ///
double JcB(double B[3], double Jc, double theta);				//Jc(B) relation Kim model
void initial_Jc(double& min);													  //initial critical current density Jc from magnetic field 
void critical_Jc();																			//critical current density Jc from magnetic field, applied damping factor l
void symm_dJ();													  							//copy dJ to the rest of sectors according symmetry
void load_Jcbdata();																		//load Jc(B) data from measurement

	/// Dissipation ///
double xdiss_UJ(double J[3], double Jc, double B[3], double N, int m);							//dissipation	
void ES_EJ(double J[3], int m, double N, double Jcx, double B[3], double E[3]);  		//electrostatic field E(J) 
double L3(int m, double J[3]);																											//calculates L3=V*(U(J)-U(J)previous)
double U(int m);																																		//elemental function, Power law
void U_dTz(int m);																																	//calculates U at edge Z
void U_dTy(int m);																																	//calculates U at edge Z
void U_dTx(int m);																																	//calculates U at edge Z
double UeZ(int m, double Jx1, double Jx2, double Jy1, double Jy2);									//calculates U according dTz
double UeX(int m, double Jy1, double Jy2, double Jz1, double Jz2);  								//calculates U according dTx
double UeY(int m, double Jx1, double Jx2, double Jz1, double Jz2);									//calculates U according dTy
void recal_U();																																			//U in cells
void recal_nB();																																		//n in surfaces according magnetic field

	///T vector ///
double Jx_Tdl(int m);																		//dJx at surface X m from 4 enclosed edges Y,Z
double Jy_Tdl(int m);																		//dJy at surface Y m from 4 enclosed edges X,Z
double Jz_Tdl(int m);																		//dJz at surface Z m from 4 enclosed edges X,Y
void cal_dJfromdT();																		//calculate dJ from dT
void symm_dT();															 						//copy dT vector to the rest of sectors according symmetry

 /// Minimization ///
void sector_iter(int nsec);															//1 itteration step of T without symmetry

	///Energy///
double dTz(int m, double ht);														//change of energy because of dTz
double dTx(int m, double ht);														//change of energy because of dTx	
double dTy(int m, double ht);														//change of energy because of dTy	
double L1_x(int m, double h);														//interaction with self-field
double L2_x(int m, double h);														//interaction with self-field	
double L1_y(int m, double h);														//interaction with self-field
double L2_y(int m, double h);														//interaction with applied-field	
double L1_z(int m, double h);														//interaction with applied-field
double L2_z(int m, double h);														//interaction with applied-field

	///Update Energy///
void update_dTx(int m, double ht, int nsec);						//updates dJ according dTx
void update_dTy(int m, double ht, int nsec);						//updates dJ according dTy
void update_dTz(int m, double ht, int nsec);						//updates dJ according dTz

	/// AC loss ///
void local_loss(int m);																	//local loss per cell
double loss();																					//total loss per period calculates between dt/T (0.25-0.75)

	/// Is current ///
void set_us();																					//set unit vector for source current us in the mesh		

	/// Other ///
void save_matrix();																			//calculates all matrix (A, H)
void save_temp();																				//save temporery variable dJ,dT to Jx, dTx
void damping(double l);																	//applied damping factor to J and average potential
void dif(double& dif, double l,int& sum);								//find difference between two iteration steps
void set_Ba(int it);																		//set instant applied field B_a, Ba
void set_Is(int it);																		//set current source Is according time line
void set_Jc(int it);																		//calculates Jc according time line
void time_variable();	 																	//Jo+=dJ, dJ=0, Ao+=dA, dA=0  

///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Save to file ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

void clear_files();																			//clear contex of all outputs files
void save_par();																				//save parameters for gnuplot
void save_3d();																					//save variables from cells n rc 0 J |J| Jcx A_a B_J B T m
void save_plane();																			//save variables from plane
void save_surface();																		//save variables to surfaces RC J  dJ			
void cal_AV();																					//calculates average of Jx,Jy,Jz and Bx,By,Bz over thickness	
void save_AV();																					//save average of Jx,Jy,Jz and Bx,By,Bz	Tx,Ty,Tz over thickness	
void save_Tz();																					//save Tz in form for current flux lines
};

///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Master class ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

class class_mcube: public class_cube 
{
	class_cube *vcube;
	int ncube;

	public:	

	void init();																					//set input variables , geometry, matrix, range 
	void cycle_dt();																			//minimization over time
};
