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
	double av_dA[3];									//delta average potential
	double av_A0[3];									//zero average potential
	double av_A[3];										//total average potential
	double dAx;							 					//temporary variable delta potential 
	double us;												//unit vector of source current
	double Is;												//source current
	double dI;												//delta of source current
	double dAs;												//delta vector potential of source current
	double As[3];											//total vector potential of source current
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
	double av_dA[3];									//delta average potential
	double av_A0[3];									//zero average potential
	double av_A[3];										//total average potential
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
	double av_dA[3];									//delta average potential
	double av_A0[3];									//zero average potential
	double av_A[3];										//total average potential
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
	double Js[3];											//total current density from current source	
	double Js0[3];										//total current density from current source in privious time step t0	
};

struct struct_sector  
{
	int adr[3];												//i,j,k address
  int *sadrY;                       //Y surface adress in the sector  
  int snY;                          //total number of Y surfaces in the sector
  int *sadrZ;                       //Z surface adress in the sector  
  int snZ;                          //total number of Z surfaces in the sector
	double dJ_av[3];			   					//average delta current density dJ_av(J_x, J_y, J_z)
	double rcx[3];			    					//vector position of X sector
	double rcy[3];			    					//vector position of Y sector
	double rcz[3];			    					//vector position of Z sector
	double Vi[3];				    					//volume of influence in sector for surfaces X,Y,Z
	double Px[3];											//dipole term for x component 
	double Py[3];											//dipole term for y component
	double Pz[3];											//dipole term for z component
	double uy;                        //unit vector
};

class class_cube      
{
	protected:

  //Input variables//	
	double x;													//size of x [m]
	double xl;												//thickness of linear material [m]
	double y;													//size of y [m]
 	double z;													//size of z [m]
	double x0;												//poistion x [m] of coordinate system 
	double y0;												//poistion y [m] of coordinate system
	double z0;												//poistion z [m] of coordinate system
	double R1;												//first radius of sample [m]
	double R2;												//second radius of sample [m] 
	double dR;												//thickness of sample [m]
	double FI;												//azimuth [degrees]
	double FI1;												//azimuth [degrees]
	double Z;													//hight [m]		
	int sys;													//coordinate system 
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
	int Jvar;													//homogenous critical current density in sample
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
	int shape;												//shape of thin film
	int num_threads;									//number of threads	
	int ijk_polar;										//i,j,k address in polar coordinate system
	int shiftx1;											//shift of sectors set2 along Z axes 			
	int shiftx2;											//shift of sectors set3 along Z axes
	int shifty1;											//shift of sectors set2 along Z axes 			
	int shifty2;											//shift of sectors set3 along Z axes
	int shiftz1;											//shift of sectors set2 along Z axes 			
	int shiftz2;											//shift of sectors set3 along Z axes
	int multi_e;											//multipole expansion of vector potential
	int nsce[3];											//number of cells in sector for evaluation of vector potential
	int ne;										  			//radius of the neigberhood sectors for vector potential evaluation 


	int nsx[3];												//numbers of sectors in X direction for three sets 
	int nsy[3];												//numbers of sectors in Y direction for three sets
	int nsz[3];												//numbers of sectors in Z direction for three sets
	int nsxe;													//numbers of sectors in X direction for av_dA evaluation
	int nsye;													//numbers of sectors in Y direction for av_dA evaluation
	int nsze;													//numbers of sectors in Z direction for av_dA evaluation
	int NCX;													//total number of cells in sectors along the X axis
	int NCY;													//total number of cells in sectors along the Y axis
	int NCZ;													//total number of cells in sectors along the Z axis	

	//sectors//
	int set;													//total number of sets set=0 not shifted, set=1 shift1, set=2 shift2	
	int nsall;												//total number of sectors (set0+set1+set2)
	int xsall[3];											//total number of not shifterd sectorsxsall[0]=set0,....
	int nseall;												//total number of sectors for av_dA evaluation
	int nscx;													//total number of X cells in one sector
	int nscy;													//total number of Y cells in one sector
	int nscz;													//total number of Z cells in one sector
  
	int *cii, *cai;										//minimum and maximum of i at cells
	int *cij, *caj;										//minimum and maximum of j at cells
	int *cik, *cak;										//minimum and maximum of k at cells
	int *ceii, *ceai;									//minimum and maximum of i at cells for multipole expansion
	int *ceij, *ceaj;									//minimum and maximum of j at cells for multipole expansion
	int *ceik, *ceak;									//minimum and maximum of k at cells for multipole expansion

	int *sxii, *sxai;									//minimum and maximum of i at cells
	int *sxij, *sxaj;									//minimum and maximum of j at cells
	int *sxik, *sxak;									//minimum and maximum of k at cells

	int *syii, *syai;									//minimum and maximum of i at cells
	int *syij, *syaj;									//minimum and maximum of j at cells
	int *syik, *syak;									//minimum and maximum of k at cells

	int *szii, *szai;									//minimum and maximum of i at cells
	int *szij, *szaj;									//minimum and maximum of j at cells
	int *szik, *szak;									//minimum and maximum of k at cells


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

	//Polar coordinate system//
	int FIdegree;											//FI in degree units
	double dr;												//dr=dR/ncx size of cell in radial direction
	double dfi;												//dfi=FI/ncy
	double dz;												//dz=Z/ncz
	double sxl;												//left most position of the sample
	double sxr;												//right most position of the sample
	double syb;												//bottom most position of the sample
	double syt;												//top most position of the sample
	double drdfi;											//difference in radial and azimuthal direction for one elementh
	int turn;													//number of turns in coil
	int nc_turn;											//number of cells in one turn

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
	struct_sector *vs;								//variabels of sectors

	double *cmAx;											//potencial matrix for surfaces X size 1 x nsurx
	double *cmAy;											//potencial matrix for surfaces Y size 1 x nsury
	double *cmAz;											//potencial matrix for surfaces Z size 1 x nsurz

	double **bcmAz1;									//potencial matrix for surfaces z size ncx x ncz
	double **bcmAz2;									//potencial matrix for surfaces z size ncx x ncz
	double **bcmAx1;									//potencial matrix for surfaces x size ncx x ncz
	double **bcmAx2;									//potencial matrix for surfaces x size ncx x ncz

	double *snAx;											//vector potencial matrix for surfaces X size 1 x nsurx
	double *snAy;											//vector potencial matrix for surfaces Y size 1 x nsury
	double *snAz;											//vector potencial matrix for surfaces Z size 1 x nsurz
	double *snAxy;										//vector potencial matrix for surfaces X size 1 x nsurxy

	double *****zmAx;									//potencial matrix for surfaces X
	double *****zmAy;									//potencial matrix for surfaces Y
	double *****zmAz;									//potencial matrix for surfaces Z
	double ****zmAxy;									//potencial matrix for surfaces Y
	double ****zmAxz;									//potencial matrix for surfaces Z

	double **mAx;											//potencial matrix for surfaces X	r
	double **mAy;											//potencial matrix for surfaces Y fi
	double **mAz;											//potencial matrix for surfaces Z z
	double **mAxy;										//potencial matrix for surfaces X r-fi

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

	double **emAx;										//vector potencial matrix for surfaces X and sectors
	double **emAy;										//vector potencial matrix for surfaces Y and sectors
	double **emAz;										//vector potencial matrix for surfaces Z and sectors

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
int ijktonedgeXp(int i,int j,int k);										//polar address of edge X n from i,j,k
int ijktonedgeY(int i,int j,int k);											//address of edge Y n from i,j,k
int ijktonedgeYp(int i,int j,int k);										//polar address of edge Y n from i,j,k
int ijktonedgeZ(int i,int j,int k);											//address of edge Z n from i,j,k
int ijktonedgeZp(int i,int j,int k);										//polar address of edge Z n from i,j,k
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
int ijktonsurXp(int i,int j,int k);											//polar address of surface X n from i,j,k
int ijktonsurYp(int i,int j,int k);											//polar address of surface Y n from i,j,k
int ijktonsurZp(int i,int j,int k);											//polar address of surface Z n from i,j,k
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
void cell_var();																				//set various critical current density in sample	
int icell(int i,int j,int k);			    		 							//address of cell n from i,j,k
int icellp(int i,int j,int k);			   		 							//polar address of cell n from i,j,k
int iplanecell(int i,int j,int k);			    		 				//address of cell n from i,j,k
void ncelltoijk(int n, int& i, int& j, int& k);         //address of cell i,j,k from n
void ncellplanetoijk(int n,int& i,int& j,int& k); 			//address of cell i,j,k from n

  ///Sectors///
void nrange(int nsec);																	//adresses of edge X,Y,Z, cell and surfaces X,Y,Z in n sector
void ne_range();																				//adresses of cells X,Y,Z in sectors for vector potential multipole expansion
void nrange_db(int nsec);																//cout adresses of cells in all sectors
void ne_range_db();																			//cout adresses of cells in all sectors
void cal_sectors();																			//set optimum cells per sector (nscx...) and calculates number of sectors (nsx....)
void cal_e_sectors();																		//calculates number of sectors for vector potential multipole expansion
void nsectortoijk(int s,int n,int& i,int& j,int& k); 		//address of sector i,j,k from n
double ijktosec(int i, int j, int k);										//address of sector n from i,j,k
void vol_sector(int n, double V[3]);										//volume of the influence for expansion sector
void sector_rc(int s, int m, double rc[3]);							//vector position of the expansion sector
void recal_e_sec();																			//recalculate dJ_av, P[9],... in the expansion sectors
void vs_sadr();																					//i,j,k address of sector for multipole expansio at surfaces
void vs_inner_box();																		//range of surfaces for inner multipole expansion
void s_dis();																						//surface distance to the sector for the multipole expansion

	/// Interoplation ///
double interX(double r[3], int n);											//interpolation of any point rc from the surface X n  
double interY(double r[3], int n);											//interpolation of any point rc from the surface Y n  
double interZ(double r[3], int n);											//interpolation of any point rc from the surface Z n 
double interR(double r[3], int n);											//interpolation of any point rc from the surface R n
double interfi(double r[3], int n);											//interpolation of any point rc from the surface fi n
double interZp(double r[3], int n);											//interpolation of any point rc from the surface Z n   
void interJ(double r[3], double J[3]);									//interpolated dJ of any point rc from surfaces   
void interJA_polar(double rc[3], double Jp[3], double Aap[3]);//interpolation of J(r,fi,z) at any point rc[r,fi,z]
void interJ_cell(int n, double J[3]);										//interpolated J of cell n (Jo+dJ) 
void interJ0_cell(int n, double J[3]);									//interpolated J0 of cell n Jo 
double interJs_cell(int n, int s);											//interpolated Js of cell n 
void intercell_A(int n, double A[3]);										//interpolated A of cell n
void intersurX_A(int n, double& A0, double& A1);				//interpolated A of surfaceX n
void intersurY_A(int n, double& A0, double& A1);				//interpolated A of surfaceY n
void intersurZ_A(int n, double& A0, double& A1);				//interpolated A of surfaceZ n
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
void save_sectors();																		//save sector's variables (adr, V, J_ab) to memory

///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Electrodynamic  ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector potential ///
double surX_dA (int i); 																//dA and A from applied field at surface X
double surY_dA (int i); 																//dA and A from applied field at surface Y
double surZ_dA (int i);							 										//dA and A from applied field at surface Z
void cal_av_dA (int s);		 															//average potential at surfaces
void multipole_expansion_av_dA();												//evaluation of vector potential by multipole expansion
void fast_matrix_cmA();																	//average vector potencial matrix of surface X, Y, Z in cartezian coordinate system
void full_matrix_mA();																	//average vector potencial matrix of surface X, Y, Z in cartezian coordinate system
void multipole_emA();																		//average vector potencial matrix of surface X, Y, Z in cartezian coordinate system
void full_polar_matrix_A();															//average vector potencial matrix of surface X, Y, Z in polar coordinate system
void fast_polar_matrix_A();															//average vector potencial matrix of surface X, Y, Z in polar coordinate system
void symZ_polar_matrix_A();															//average vector potencial matrix of surface X, Y, Z in polar coordinate system
void full_polar_matrix_A_by_symZ();											//average vector potencial matrix of surface X, Y, Z in polar coordinate system, matrix fill by zmAx...
double read_cmAx(int j, int i);													//read from cmAx (vector potential matrix)		
double read_cmAy(int j, int i);													//read from cmAy (vector potential matrix)
double read_cmAz(int j, int i);													//read from cmAz (vector potential matrix)
double read_cmAxy(int i, int j);												//read from cmAxy (vector potential matrix)	
double read_cmAxz(int i, int j);												//read from cmAxy (vector potential matrix)		
double num_a_auto(int l, int m, int i); 								//automaticly find the best resolution of sub-elemets
double num_a(int l, int m, int na[3], int nb[3],int s);	//average scalar potential of volume created by volume of 1 element, numerical calculation
double num_ar(int l,int m, int na[3], int nb[3],int s);	//scalar potential of volume created 1 element, numerical calculation
double num_arfi(int l,int m, int na[3], int nb[3],int s);//scalar potential of volume created 1 element, numerical calculation
double num_arfi1(int l,int m, int na[3], int nb[3],int s);//scalar potential of volume created 1 element, numerical calculation
double num_arz(int l,int m, int na[3], int nb[3],int s);//scalar potential of volume created 1 element, numerical calculation
void boundary_cmA();																		//average vector potencial matrix of surface X,Y in cartezian coordinate system
void num_Ay_boundary(int& i, int& j);										//shits adresses of the surfaces Y in boundary from second end to the first one (issue with condition for creating overlapping mesh - B,Al)	
int num_Ay_boundary1(int i);														//shits adresses of the surfaces Y in boundary from second end to the first one (issue with condition for creating overlapping mesh - B,Al)	
void matrix_snA();																			//average vector potencial for source current

	/// Magnetic field ///
void magnetic_field();																	//recalculate magnetic field B= B_J(J) + Ba
void external_magnetic_field();													//recalculate external magnetic field B= B_J(J) + Ba
void magnetic_moment(int n, double m[3]);								//magnetic moment at the cell
double fast_matrix_Hx();																//average magnetic field matrices cmHx, cmHy, cmHz
double full_matrix_Hx();																//average magnetic field matrices cmHx, cmHy, cmHz
void symZ_polar_matrix_H();															//magnetic field matrix of surface X, Y, Z in polar coordinate system
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
double num_Hr(int m, int l, int na[3], int nb[3],int s, int t);//magnetic field of volume created by volume of 1 element, numerical calculation


 /// Current density ///
double J_flux();																				//flux lines of current
double JcB(double B[3], double Jc, double J[3], double theta);//Jc(B) relation Kim model
void initial_Jc(double& min);													  //initial critical current density Jc from magnetic field 
void critical_Jc();																			//critical current density Jc from magnetic field, applied damping factor l
void critical_Jc1();
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
void dif(double& dif, double l,int& sx,int& sy,int& sz);//find difference between two iteration steps
void set_Ba(int it);																		//set instant applied field B_a, Ba
void set_Is(int it);																		//set current source Is according time line
void set_Jc(int it);																		//calculates Jc according time line
void time_variable();	 																	//Jo+=dJ, dJ=0, Ao+=dA, dA=0 
void cartezian_polar(double c[3],double p[3]);					//transform cartezian coorinates to polar
void polar_cartezian(double p[3],double c[3]);					//transform polar coorinates to cartezain  

///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Save to file ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

void clear_files();																			//clear contex of all outputs files
void save_par();																				//save parameters for gnuplot
void save_3d();																					//save variables from cells n rc 0 J |J| Jcx A_a B_J B T m
void save_3d_cross();																		//save variables from cells to file format for multiplanes format RC[3] J[3] B[3]
void save_plane();																			//save variables from plane
void save_surface();																		//save variables to surfaces RC J  dJ
void save_mHx();																				//save matrices mHzxy....
void save_cmAx();																				//save matrices cmAx....
void load();																						//load output3Dz.txt
void load1();																						//load output3Dy.txt
void save_edge();																				//save variables from edges n rc dt w wi				
void cal_AV();																					//calculates average of Jx,Jy,Jz and Bx,By,Bz over thickness	
void save_AV();																					//save average of Jx,Jy,Jz and Bx,By,Bz	Tx,Ty,Tz over thickness	
void save_Tz();																					//save Tz in form for current flux lines

void recal_polar_data();																//transform data from polar coordiante system to cartezian coordinate system

void test();																						//Ax
void test1();																						//axij
void check();
void load_matrix();																			//load matrix cmAx/y/xy
void load_matrix1();																		//load matrix cmAx/y/xy
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



