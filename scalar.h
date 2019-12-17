using namespace std;

////////////////////////////////////////////
//////////////////point/////////////////////
////////////////////////////////////////////

//scalar potential of point 
double spp(double *r, double *r0, double q)		
{
	double pi;
	double ep=8.854187817620e-12;
	double rpp[3];
	pi=4.0*atan(1.0);
	vec_subs2(r,r0,rpp);
	return (q/(4.0 * pi * ep * vec_mod(rpp) ));			
}

//vector electric field of point
void vefp (double *r, double *r0, double q, double *E)
{	
	double pi;
	double ep=8.854187817620e-12;
	double rpp[3];
	int i;
	pi=4.0*atan(1.0);
	vec_subs2(r,r0,rpp);
		for (i=0;i<=2;i++) 
				E[i]=((q * (rpp[i])) / (4 * pi * ep * pow(vec_mod(rpp),3.0) ));
}

////////////////////////////////////////////
//////////////////line//////////////////////
////////////////////////////////////////////

//scalar potential of line 
double spl (double *r, double *r0, double q, double l)							
{
	double pi, x1, x2;
	double ep=8.854187817620e-12;
	double rpp[3];
	pi=4.0*atan(1.0);
	vec_subs2(r,r0,rpp);
  x1= rpp[0]+(l/2);
	x2= rpp[0]-(l/2);	

  if ( (rpp[1]==0) && (rpp[2]==0) )  return ((q/l)/(4.0*pi*ep))*(((x1/fabs(x1))*log(2*fabs(x1))) - (x2/fabs(x2))*log(2*fabs(x2)));
	if (x1==l/2) return ((q/l)/(4.0*pi*ep))*((   asinh((x2)/(sqrt(rpp[1]*rpp[1] + rpp[2]*rpp[2])))));
	if (x2==-l/2)  return ((q/l)/(4.0*pi*ep))*(( - asinh((x1)/(sqrt(rpp[1]*rpp[1] + rpp[2]*rpp[2])))));	
	
	return ((q/l)/(4.0*pi*ep))*((asinh((x1)/(sqrt(rpp[1]*rpp[1] + rpp[2]*rpp[2])))) - (asinh((x2)/(sqrt(rpp[1]*rpp[1] + rpp[2]*rpp[2]))))); 

	if (rpp[0]==0) return 0;
}

//vector electric field of line
void vefl (double *r, double *r0, double q, double l, double *E)
{	
	double pi,a,c,x1,x2,y,z;
	double ep=8.854187817620e-12;
	double rpp[3];
	pi=4.0*atan(1.0);
	vec_subs2(r,r0,rpp);

	c= (q/l)/(4*pi*ep);
	x1= rpp[0]+(l/2);
	x2= rpp[0]-(l/2);
	a = (pow(rpp[1],2) + pow(rpp[2],2));
	E[0]= (-1) * c * ( (1/sqrt(a + x1*x1)) - (1/sqrt(a + x2*x2)) );
  E[1]= (-1) * c * ( ((x2*rpp[1])/((sqrt(a + x2*x2))*a)) - ((x1*rpp[1]) / ((sqrt(a + x1*x1))*a))  );
	E[2]= (-1) * c * ( ((x2*rpp[2])/((sqrt(a + x2*x2))*a)) - ((x1*rpp[2]) / ((sqrt(a + x1*x1))*a))  );
}

////////////////////////////////////////////
/////////////////surface////////////////////
////////////////////////////////////////////

//scalar potential of surface 
double sps(double *r, double *r0, double q, double a, double b)			
{
	double pi, c, z;
	double ep=8.854187817620e-12;
	double rpp[3];
	double phi=0, phi1=0;
	double mFS[2][2]={0}, I[2][2]={0};
	double FS (double xpp,double ypp, double rpp);
	int i,j;	
	vec_subs2(r,r0,rpp);
	pi=4.0*atan(1.0);
	c  = (q/(a*b))/(4*pi*ep);		
	z  = rpp[2];
  I[0][0] = rpp [0]+(a/2); 	//x1
	I[0][1] = rpp [0]-(a/2);	//x2
	I[1][0] = rpp [1]+(b/2);	//y1
	I[1][1] = rpp [1]-(b/2);	//y2

//loop to evaluate fundamental function FS with boundary
	for(i=0; i<=1; i++)			
   { 
			for (j=0; j<=1; j++)	
		    {			mFS[i][j] = FS ( I[0][i] , I[1][j], z);
			   			mFS[i][j] = pow(-1,i+j)*mFS[i][j];
		          phi1+=mFS[i][j]; 
				}
	 }	
	phi= c * phi1;
	return phi;	
}

//fundamental function scalar potential of surface without boundary
double FS (double xpp, double ypp, double z)
{

	if ((xpp==0)&&(z==0)||(ypp==0)&&(z==0)) return 0;
	else	

 	return /*1*/ypp*asinh(xpp/sqrt(ypp*ypp + z*z)) - /*2*/z*atan(xpp*ypp/(z*sqrt(xpp*xpp + ypp*ypp + z*z))) + /*3*/xpp*asinh(ypp/(sqrt(xpp*xpp + z*z)));
}

//vector electric field of surface 
void vefs (double *r, double *r0, double q, double a, double b, double *E)
{	
	double FES (int i, double xpp, double ypp, double z);	
	double pi,c,z,x1,x2,y1,y2;
	double ep=8.854187817620e-12;
	double rpp[3]={0};
	int i;
	pi=4.0*atan(1.0);
	vec_subs2(r,r0,rpp);

	c  = (q/a*b)/(4*pi*ep);
	z  = rpp[2];
	x1 = rpp [0]+(a/2); 
	x2 = rpp [0]-(a/2);	
	y1 = rpp [1]+(b/2);	
	y2 = rpp [1]-(b/2);

for(i=0 ; i<=2 ; i++)
	{
	  E[i]=(-1) * c * ( FES (i, x1, y1, z) - FES (i, x1, y2, z) - FES (i, x2, y1, z) + (FES (i, x2, y2, z)) );  
	}
}

//fundamental function vector electric field of surface without boundary
double FES (int n, double xpp, double ypp, double z)
{
	double phi,x2,y2,z2,a,b;
	x2 = xpp*xpp;
	y2 = ypp*ypp;
	z2 = z*z;
	a  = x2 + y2 + z2;
	b  = sqrt (a);

switch (n) {
  case 0:
			if ((xpp==0)&&(ypp==0)&&(z==0)) return 0;
			if ((xpp==0)&&(ypp==0)) return (/*3*/log(ypp + b));
			if ((xpp==0)&&(z==0)) return (/*1*/ypp/fabs(ypp) + /*3*/(ypp/fabs(ypp))*log(2*(fabs(ypp))));
			if ((ypp==0)&&(z==0)) return (/*3*/log(ypp + b));
			if (xpp==0) return (/*1*/(ypp/b) - /*2*/(ypp/(sqrt(y2+z2))) + /*3*/log(ypp + b));
		  if (ypp==0) return (/*3*/log(ypp + b));
			if (z==0) return (/*1*/(ypp/b) + /*3*/log(ypp + b) - /*4*/(ypp/sqrt(x2 + y2))); 

			return /*1*/ypp/b - /*2*/ypp*z2/((x2+z2)*b) + /*3*/log(ypp + b) - /*4*/x2*ypp/(b*(x2 + z2));
    break;

  case 1:
			if ((xpp==0)&&(ypp==0)&&(z==0)) return 0;
			if ((xpp==0)&&(ypp==0)) return /*1*/log(xpp + b);
			if ((xpp==0)&&(z==0)) return /*1*/log(xpp + b);
			if ((ypp==0)&&(z==0)) return (/*1*/(xpp/fabs(xpp))*log(2*(fabs(xpp))) + /*4*/(xpp/fabs(xpp)));
      if (xpp==0) return /*1*/log(xpp + b);
			if (ypp==0) return /*1*/log(xpp + b);
  	  if (z==0) return /*1*/log(xpp + b);

    	return /*1*/log(xpp + b) - /*2*/xpp*y2/(b*(y2 + z2)) - /*3*/xpp*z2/((y2+z2)*b) + /*4*/(xpp/b); 
    break;

  case 2: 
			if ((xpp==0)&&(ypp==0)&&(z==0)) return 0;
			if ((xpp==0)&&(ypp==0)) return 0;
			if ((xpp==0)&&(z==0)) return 0;
			if ((ypp==0)&&(z==0)) return 0;
			if ((xpp==0)||(ypp==0)||(z==0)) return 0;

			return - /*1*/((xpp*ypp*z)/(b*(y2 + z2))) - /*2*/atan((xpp*ypp)/(z*b)) + /*3*/((xpp*ypp*z)*(x2 + y2 + 2*z2)/((x2+z2)*(y2+z2)*b)) - /*4*/((xpp*ypp*z)/(b*(x2 + z2))); 
	break;
					}
}

////////////////////////////////////////////
//////////////////volume////////////////////
////////////////////////////////////////////

//scalar potential of volume
double spv(double *r, double *r0, double q, double a, double b, double c)			
{
  double FV (double xpp, double ypp, double zpp);
	double pi,e,x1,x2,y1,y2,z1,z2;
	double ep=8.854187817620e-12, mi0;
	double rpp[3]={0};
	double phi=0, phi1=0;
	double mFV[2][2][2]={0}, I[3][2]={0};
	int i,j,k;			
	pi=4.0*atan(1.0);
	vec_subs2(r,r0,rpp);

	mi0=4.0*pi*1e-7;

	e  = mi0/(a*b*c*4*pi);

//	e  = (q/(a*b*c))/(4*pi*ep);		

 	I[0][0] = rpp [0]+(a/2); //x1
	I[0][1] = rpp [0]-(a/2); //x2	
	I[1][0] = rpp [1]+(b/2); //y1	
	I[1][1] = rpp [1]-(b/2); //y2	
	I[2][0] = rpp [2]+(c/2); //z1
	I[2][1] = rpp [2]-(c/2); //z2

//loop to evaluate fundamental function FV with boundary
for(i=0; i<=1; i++)
   {for (j=0; j<=1; j++)
				{for (k=0; k<=1; k++)
			       {		mFV[i][j][k] = FV ( I[0][i], I[1][j], I[2][k] );
									mFV[i][j][k] = pow(-1,i+j+k)*mFV[i][j][k];	  
									phi1+= mFV[i][j][k];
						 }
 			  }	
	 }
	phi= e * phi1;
	return phi;
}

//fundamental function scalar potential of volume without boundary
double FV (double xpp, double ypp, double zpp)
{
	double a,b,c,phi,x2,y2,z2;

	x2 = xpp*xpp;
	y2 = ypp*ypp;
	z2 = zpp*zpp;
	a = sqrt(y2 + z2);
	b = sqrt(x2 + z2);
 	c = sqrt(x2 + y2 + z2);

  if ((xpp==0)&&(ypp==0)&&(zpp==0)) return 0;
	if ((xpp==0)&&(ypp==0)) return 0; 
	if ((xpp==0)&&(zpp==0)) return 0;
  if ((zpp==0)&&(ypp==0)) return 0; 
  if (xpp==0) return /*2*/ypp*zpp*log(xpp + c);
	if (ypp==0) return /*1*/xpp*zpp*log(ypp +c) + /*4*/(z2/2)*atan((ypp*c + y2 + z2)/(xpp*zpp));	
	if (zpp==0) return + /*3*/xpp*ypp*log(zpp + c) + /*5*/(y2/2)*atan((zpp*c + y2 + z2)/(xpp*ypp));	
 	
 	return  /*1*/xpp*zpp*log(ypp +c) + /*2*/ypp*zpp*log(xpp + c) + /*3*/xpp*ypp*log(zpp + c) + /*4*/(z2/2)*atan((ypp*c + y2 + z2)/(xpp*zpp)) + /*5*/(y2/2)*atan((zpp*c + y2 + z2)/(xpp*ypp)) 
        - /*6*/z2*atan(xpp*ypp/(zpp*c)) - /*7*/y2*atan(xpp*zpp/(ypp*c)) - /*8*/(x2/2)*atan(ypp*zpp/(xpp*c));
}

//vector electric field of volume

void vefv (double *r, double *r0, double q, double a, double b, double c, double *E)
{	
	double FEV(int i, double xpp, double ypp, double zpp);
	double pi,d,x1,x2,y1,y2,z1,z2;
	double ep=8.854187817620e-12;
	double rpp[3]={0};
	int i;
	pi=4.0*atan(1.0);
	vec_subs2(r,r0,rpp);

	d  = (q/a*b*c)/(4*pi*ep);
	x1 = rpp [0]+(a/2); 
	x2 = rpp [0]-(a/2);	
	y1 = rpp [1]+(b/2);	
	y2 = rpp [1]-(b/2);
	z1 = rpp [2]+(c/2);
	z2 = rpp [2]-(c/2);		

//loop to evaluate fundamental function FEV with boundary

for(i=0 ; i<=2 ; i++)		

	{
		E[i]=(-1) * d * ( FEV(i, x1, y1, z1) - FEV(i, x1, y1, z2) - FEV(i, x1, y2, z1) + FEV(i, x1, y2, z2) - FEV(i, x2, y1, z1) + FEV(i, x2, y1, z2) + FEV(i, x2, y2, z1) - FEV(i, x2, y2, z2));
	}
}	

//fundamental function vector electric field of volume without boundary
double FEV (int n, double xpp, double ypp, double zpp)
{
	double phi,x2,y2,z2, x3, y3, z3, b;
	x2 = pow(xpp,2.0); x3 = pow(xpp,3.0);
	y2 = pow(ypp,2.0); y3 = pow(ypp,3.0);
	z2 = pow(zpp,2.0); z3 = pow(zpp,3.0);
	b  = sqrt (x2 + y2 + z2);

switch (n) {
  case 0:
			if ((xpp==0)&&(ypp==0)&&(zpp==0)) return 0;
			if ((xpp==0)&&(ypp==0)) return /*1*/(zpp*log(ypp + b));
			if ((xpp==0)&&(zpp==0)) return /*4*/(ypp*log(zpp + b));
			if ((ypp==0)&&(zpp==0)) return 0;
	    if (xpp==0) return /*1*/(zpp*log(ypp + b)) + /*3*/(ypp*zpp/b) + /*4*/(ypp*log(zpp + b)) - /*6*/(ypp*z2*zpp/(2*b*(x2 + z2))) - /*7*/(y2*ypp*zpp/(2*b*(x2 + y2)));
		  if (ypp==0) return /*1*/(zpp*log(ypp + b));
			if (zpp==0) return /*4*/(ypp*log(zpp + b));

			return /*1*/(zpp*log(ypp + b)) - /*2*/((x2*ypp*zpp)/(b*(x2 + z2))) + /*3*/(ypp*zpp/b) + /*4*/(ypp*log(zpp + b)) - /*5*/(x2*ypp*zpp/(b*(x2 + y2)))
  		     - /*6*/(ypp*z3/(2*b*(x2 + z2))) - /*7*/(y3*zpp/(2*b*(x2 + y2))) - /*8*/(xpp*atan(ypp*zpp/(xpp*b))) + /*9*/(x2*ypp*zpp*((2*x2 + y2 + z2)/(2*(x2 + y2)*(x2 + z2)*b)));
    break;
  case 1:
			if ((xpp==0)&&(ypp==0)&&(zpp==0)) return 0;
			if ((xpp==0)&&(ypp==0)) return /*2*/zpp*log(xpp + b);
			if ((xpp==0)&&(zpp==0)) return 0;
			if ((ypp==0)&&(zpp==0)) return /*4*/xpp*log(zpp + b);
	    if (xpp==0) return /*2*/zpp*log(xpp + b);
		  if (ypp==0) return /*1*/(xpp*zpp/b) + /*2*/zpp*log(xpp + b) +/*4*/xpp*log(zpp + b) - /*6*/(xpp*z3/(2*b*(y2 + z2))) - /*9*/x3*zpp/(2*b*(x2 + y2));
			if (zpp==0) return /*4*/xpp*log(zpp + b);

    	return  /*1*/(xpp*zpp/b) + /*2*/zpp*log(xpp + b) - /*3*/(xpp*y2*zpp/(b*(y2 + z2))) + /*4*/xpp*log(zpp + b) - /*5*/(xpp*y2*zpp/(b*(x2 + y2))) - /*6*/(xpp*z3/(2*b*(y2 + z2)))
     		    - /*7*/ypp*atan(xpp*zpp/(ypp*b)) + /*8*/xpp*y2*zpp*((x2 + 2*y2 + z2)/(2*(x2 + y2)*(y2 + z2)*b)) - /*9*/x3*zpp/(2*b*(x2 + y2));     
	  break;
  case 2: 
    	if ((xpp==0)&&(ypp==0)&&(zpp==0)) return 0;
			if ((xpp==0)&&(ypp==0)) return 0;
			if ((xpp==0)&&(zpp==0)) return /*3*/ypp*log(xpp + b);
			if ((ypp==0)&&(zpp==0)) return /*1*/xpp*log(ypp + b);
	    if (xpp==0) return /*3*/ypp*log(xpp + b);
		  if (ypp==0) return /*1*/xpp*log(ypp + b);
			if (zpp==0) return /*1*/xpp*log(ypp + b) + /*3*/ypp*log(xpp + b) + /*5*/(xpp*ypp/b) - /*8+9*/xpp*ypp/sqrt(x2+y2);
	
			return /*1*/xpp*log(ypp + b) - /*2*/(xpp*ypp*z2/(b*(x2 + z2))) + /*3*/ypp*log(xpp + b) - /*4*/(xpp*ypp*z2/(b*(y2 + z2))) + /*5*/(xpp*ypp/b) 
	         - /*6*/(zpp*atan(xpp*ypp/(zpp*b))) + /*7*/(xpp*ypp*z2*((x2 + y2 + 2*z2)/(2*(x2 + z2)*(y2 + z2)*b))) - /*8*/(xpp*y3/(2*b*(y2 + z2))) - /*9*/(x3*ypp/(2*b*(x2 + z2)));
	break;
					}
}

//vector potential of volume

double vpv (double *r, double *rs, double a, double b, double c, int s)
{	
	double FVP(int i, double xpp, double ypp, double zpp);
	double FV(double xpp, double ypp, double zpp);
	double mi0, pi,d, a1, b1;
	double rpp[3]={0}, I[3][2]={0}, mFV[2][2][2]={0};
	double rcm[3],rcp[3];
	int i=0,j=0,k=0;
	pi=4.0*atan(1.0);
	mi0=4.0*pi*1e-7;
	
	rcm[0]=rs[0]-a/2.0;
	rcm[1]=rs[1];
	rcm[2]=rs[2];

	rcp[0]=rs[0]+a/2.0;
	rcp[1]=rs[1];
	rcp[2]=rs[2];

	d  = mi0/(a*b*c*4*pi);

//loop to evaluate fundamental function FVP with boundary

	vec_subs2(r,rcm,rpp);
 	I[0][0] = rpp [0]+(a/2); //x1
	I[0][1] = rpp [0]-(a/2); //x2	
	I[1][0] = rpp [1]+(b/2); //y1	
	I[1][1] = rpp [1]-(b/2); //y2	
	I[2][0] = rpp [2]+(c/2); //z1
	I[2][1] = rpp [2]-(c/2); //z2

//loop to evaluate fundamental function FV with boundary
//left triangel
for(i=0; i<=1; i++)
	 {
		  for(j=0; j<=1; j++)
				{
					for(k=0; k<=1; k++)
						 {		
									mFV[i][j][k] = ((r[0] - (rcm[0] - (a/2)))/a) * FV ( I[0][i], I[1][j], I[2][k] ) - (1/a) * FVP (s, I[0][i], I[1][j], I[2][k] );
									mFV[i][j][k] = pow(-1,i+j+k) * mFV[i][j][k];	  
									a1+= mFV[i][j][k];
						 }
 			  }	
	 }

	vec_subs2(r,rcp,rpp);
 	I[0][0] = rpp [0]+(a/2); //x1
	I[0][1] = rpp [0]-(a/2); //x2	
	I[1][0] = rpp [1]+(b/2); //y1	
	I[1][1] = rpp [1]-(b/2); //y2	
	I[2][0] = rpp [2]+(c/2); //z1
	I[2][1] = rpp [2]-(c/2); //z2

//right triangel
for(i=0; i<=1; i++)
	 {
		  for(j=0; j<=1; j++)
				{
				  for(k=0; k<=1; k++)
						 {	
									mFV[i][j][k] = ((-r[0]+rcp[0]+(a/2) )/a) * FV ( I[0][i], I[1][j], I[2][k] ) + (1/a) * FVP (s, I[0][i], I[1][j], I[2][k] );
									mFV[i][j][k] = pow(-1,i+j+k) * mFV[i][j][k];	  
									b1+= mFV[i][j][k];
						 }
 			  }	
	 }

return d*(a1 + b1);
}	

//fundamental function vector potential of volume without boundary
double FVP (int m, double xpp, double ypp, double zpp)
{
	double phi,x2,y2,z2, x3, y3, z3, a, ax, ay, az, n;
	x2 = pow(xpp,2.0); x3 = pow(xpp,3.0);
	y2 = pow(ypp,2.0); y3 = pow(ypp,3.0);
	z2 = pow(zpp,2.0); z3 = pow(zpp,3.0);
	a  = sqrt (x2 + y2 + z2);
	ax=fabs(xpp);ay=fabs(ypp);az=fabs(zpp);
	n=1e-16; //range of zero point (okolie bodu)

switch (m) {
  case 0:

    	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
			if ((ax<n)&&(ay<n)) return /*5*/(z3/6)*log(ypp + a);
			if ((ax<n)&&(az<n)) return /*2*/(ypp/4)*(x2 + y2)*log(zpp + a);
			if ((ay<n)&&(az<n)) return 0;
	    if (ax<n) return /*1*/(ypp*zpp*a/3) + /*2*/(ypp/4)*(x2 + y2)*log(zpp + a) + /*5*/(z3/6)*log(ypp + a) - /*7*/(ypp/12)*(3*x2 + y2)*log(zpp + a);
		  if (ay<n) return /*4*/(x2*zpp/2)*log(ypp + a) + /*5*/(z3/6)*log(ypp + a);
			if (az<n) return /*2*/(ypp/4)*(x2 + y2)*log(zpp + a) + /*6*/(x2*ypp/2)*log(zpp + a) - /*7*/(ypp/12)*(3*x2 + y2)*log(zpp + a);

//			return /*1*/(ypp*zpp*a/3) + /*2*/(ypp/4)*(x2 + y2)*log(zpp + a) - /*3*/(x3/3)*atan(ypp*zpp/(xpp*a)) + /*4*/(x2*zpp/2)*log(ypp + a) + /*5*/(z3/6)*log(ypp + a) + /*6*/(x2*ypp/2)*log(zpp + a) 
//					 - /*7*/(ypp/12)*(3*x2 + y2)*log(zpp + a);

			return /*1*/(ypp*zpp*a/3) + /*2*/(y3/6)*log(zpp + a) - /*3*/(x3/3)*atan(ypp*zpp/(xpp*a)) + /*4*/(x2*zpp/2)*log(ypp + a) + /*5*/(z3/6)*log(ypp + a) + /*6*/(x2*ypp/2)*log(zpp + a);

 /*   break;
  case 1:
			if ((xpp==0)&&(ypp==0)&&(zpp==0)) return 0;
			if ((xpp==0)&&(ypp==0)) return ;
			if ((xpp==0)&&(zpp==0)) return ;
			if ((ypp==0)&&(zpp==0)) return ;
	    if (xpp==0) return ;
		  if (ypp==0) return ;
			if (zpp==0) return ;

    	return;     
	  break;
  case 2: 
    	if ((xpp==0)&&(ypp==0)&&(zpp==0)) return 0;
			if ((xpp==0)&&(ypp==0)) return 0;
			if ((xpp==0)&&(zpp==0)) return ;
			if ((ypp==0)&&(zpp==0)) return ;
	    if (xpp==0) return ;
		  if (ypp==0) return ;
			if (zpp==0) return ;
	
			return ;
	break;
*/
					}
}

/////////////////////////////////////////////////////////
//////////////////surface by surface ////////////////////
/////////////////////////////////////////////////////////

//scalar potential of surface created by surface
double spss(double *r, double *r0, double q, double a, double b,  double a1, double b1)			
{
	double pi, c, z, m=0,n=0;
	double ep=8.854187817620e-12;
	double rpp[3]={0};
	double phi=0, phi1=0;
	double FSS (double xpp,double ypp, double z);
	double mFSS[4][4]={0}, I[2][4]={0};	
	int i,j;
	pi=4.0*atan(1.0);
	vec_subs2(r,r0,rpp);

	c  = (q/a*b*a1*b1)/(4*pi*ep);		
	z  = rpp[2];
  I[0][0] = rpp [0] + (a1/2) + (a/2); //x1
	I[0][1] = rpp [0] + (a1/2) - (a/2);	//x2
	I[0][2] = rpp [0] - (a1/2) + (a/2);	//x3
	I[0][3] = rpp [0] - (a1/2) - (a/2); //x4
	I[1][0] = rpp [1] + (b1/2) + (b/2);	//y1
	I[1][1] = rpp [1] + (b1/2) - (b/2); //y2
	I[1][2] = rpp [1] - (b1/2) + (b/2);	//y3
	I[1][3] = rpp [1] - (b1/2) - (b/2); //y4

//loop to evaluate fundamental function FSS with boundary

for(i=0; i<=3; i++)
   { 
			for (j=0; j<=3; j++)
				{			mFSS[i][j] = FSS ( I[0][i], I[1][j], z);	
							if ((i == 0)||(i == 3)){m=0;}
									else m=1;					
							if ((j == 0)||(j == 3)){n=0;}
									else n=1;
  						mFSS[i][j] = pow(-1,(m+n))*mFSS[i][j];							
							phi1+=mFSS[i][j];
				}
 		}	
	phi = c * phi1;
	return phi;
}

//fundamental function scalar potential of surface by surface without boundary
double FSS (double xpp, double ypp, double zpp)
{
	double phi1, x2, y2, z2, a, a3, b, c;
	x2 = xpp*xpp;
	y2 = ypp*ypp;
	z2 = zpp*zpp;
	a  = sqrt(x2 + y2 + z2);
	a3 = pow(a,3.0);
	b  = sqrt(y2 + z2);
	c  = sqrt(x2 + z2);

    if ((xpp==0)&&(ypp==0)&&(zpp==0)) return 0;
  	if ((xpp==0)&&(ypp==0)) return (+ /*4*/((z2/2)*a) - /*6*/(a3/6)); 
  	if ((xpp==0)&&(zpp==0)) return (- /*6*/((a3)/6));
    if ((ypp==0)&&(zpp==0)) return (- /*6*/((a3)/6)); 
    if (xpp==0) return + /*4*/((z2/2)*a) - /*5*/((ypp*z2)/4)*(log((ypp + a)/(-ypp + a))) - /*6*/(a3/6);
  	if (ypp==0) return - /*1*/((xpp*z2)/4)*(log((xpp + a)/(-xpp + a))) + /*4*/((z2/2)*a) - /*6*/(a3/6);	
    if (zpp==0) return + /*2*/(xpp*y2/2)*log(xpp + a) - /*6*/(a3/6) + /*7*/(x2*ypp/2)*log(ypp + a);	

	return - /*1*/((xpp*z2)/4)*(log((xpp + a)/(-xpp + a))) + /*2*/(xpp*y2/2)*log(xpp + a) - /*3*/(xpp*ypp*zpp)*atan((xpp*ypp)/(zpp*a)) 
	       + /*4*/((z2/2)*a) - /*5*/((ypp*z2)/4)*(log((ypp + a)/(-ypp + a))) - /*6*/(a3/6) + /*7*/(x2*ypp/2)*log(ypp + a);	
}

//vector electric field of surface by surface, q total charge 
double vefss (double *r, double *r0, double q, double a, double b, double a1, double b1, int l)
{	
	double FESS(int l, double xpp, double ypp, double z);
	double pi,phi=0, d, z;
	double ep=8.854187817620e-12;
	double rpp[3]={0},I[2][4]={0},mFESS[4][4]={0} ;
	int i, j, k, m, n, o;
	pi=4.0*atan(1.0);
	vec_subs2(r,r0,rpp);

	d  = (q/(a*b*a1*b1))/(4*pi*ep);
	z  = rpp[2];
  I[0][0] = rpp [0] + (a1/2) + (a/2); //x1
	I[0][1] = rpp [0] + (a1/2) - (a/2);	//x2
	I[0][2] = rpp [0] - (a1/2) + (a/2);	//x3
	I[0][3] = rpp [0] - (a1/2) - (a/2); //x4
	I[1][0] = rpp [1] + (b1/2) + (b/2);	//y1
	I[1][1] = rpp [1] + (b1/2) - (b/2); //y2
	I[1][2] = rpp [1] - (b1/2) + (b/2);	//y3
	I[1][3] = rpp [1] - (b1/2) - (b/2); //y4

//loop to evaluate fundamental function FESS with boundary

phi=0;
	for(i=0; i<=3; i++)
 		{ 
			for (j=0; j<=3; j++)
			  {		 		
								mFESS[i][j] = FESS (l, I[0][i], I[1][j], z);	
		         		if ((i == 0)||(i == 3)){m=0;}
			         		else m=1;					
				      	if ((j == 0)||(j == 3)){n=0;}
					        else n=1;
								mFESS[i][j] = pow(-1,(m+n))*mFESS[i][j];
			          phi+=mFESS[i][j];
			  }
		}	
  return d * (-1) * phi; 
}

//magnetic field of surface by surface 
double magfss (double *r, double *r0, double a, double b, double a1, double b1, int l)
{	
	double FESS(int l, double xpp, double ypp, double z);
	double mi0, pi,phi=0, d, z;
	double ep=8.854187817620e-12;
	double rpp[3]={0},I[2][4]={0},mFESS[4][4]={0} ;
	int i, j, k, m, n, o;
	pi=4.0*atan(1.0);
	mi0=4.0*pi*1e-7;
	vec_subs2(r,r0,rpp);

	d  = (mi0/(a*b*a1*b1))/(4*pi);
	
	z  = rpp[2];

  I[0][0] = rpp [0] + (a1/2) + (a/2); //x1
	I[0][1] = rpp [0] + (a1/2) - (a/2);	//x2
	I[0][2] = rpp [0] - (a1/2) + (a/2);	//x3
	I[0][3] = rpp [0] - (a1/2) - (a/2); //x4
	I[1][0] = rpp [1] + (b1/2) + (b/2);	//y1
	I[1][1] = rpp [1] + (b1/2) - (b/2); //y2
	I[1][2] = rpp [1] - (b1/2) + (b/2);	//y3
	I[1][3] = rpp [1] - (b1/2) - (b/2); //y4

//loop to evaluate fundamental function FESS with boundary

phi=0;
	for(i=0; i<=3; i++)
 		{ 
			for (j=0; j<=3; j++)
			  {		 		
								mFESS[i][j] = FESS (l, I[0][i], I[1][j], z);	
		         		if ((i == 0)||(i == 3)){m=0;}
			         		else m=1;					
				      	if ((j == 0)||(j == 3)){n=0;}
					        else n=1;
								mFESS[i][j] = pow(-1,(m+n))*mFESS[i][j];
			          phi+=mFESS[i][j];
			  }
		}	
  return d * (-1) * phi; 
}
			
//fundamental function vector electric field of surface by surface without boundary
double FESS (int l, double xpp, double ypp, double z)
{
	double x2, y2, z2, x3, y3, z3, a, a2, a3, phi, ax, ay, az, n;
	x2 = pow(xpp,2.0);      y2 = pow(ypp,2.0);      z2 = pow(z,2.0);	
	x3 = pow(xpp,3.0);      y3 = pow(ypp,3.0);      z3 = pow(z,3.0);
	a  = sqrt (x2 + y2 + z2); a2 = pow(a,2.0);      a3 = pow(a,3.0);

	ax=fabs(xpp);ay=fabs(ypp);az=fabs(z);
	n=1e-16; //range of zero point (okolie bodu)

switch (l) {
  case 0:

  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
   	if ((ax<n)&&(ay<n)) return /*1*/ -(z2/4)*log((xpp + a)/(-xpp + a));
  	if ((ax<n)&&(az<n)) return + /*3*/ (y2/2)*log((xpp + a)/(sqrt(y2 + z2)));
    if ((ay<n)&&(az<n)) return - /*9*/ (xpp*a/2);
    if (ax<n) return - /*1*/ (z2/4)*log((xpp + a)/(-xpp + a)) + /*3*/ (y2/2)*log((xpp + a)/(sqrt(y2 + z2)));
   	if (ay<n) return - /*1*/ (z2/4)*log((xpp + a)/(-xpp + a)) - /*2*/ (xpp*z2)/(2*a) + /*7*/ ((z2*xpp)/(2*a)) - /*9*/ (xpp*a/2);
    if (az<n) return + /*3*/ (y2/2)*log((xpp + a)/(sqrt(y2 + z2))) + /*4*/ xpp*y2/(2*a) - /*9*/ (xpp*a/2) + /*10*/ (xpp*ypp)*log((ypp + a)/(sqrt(x2+z2))) 
										 - /*11*/ (x3*y2/(2*a*(x2+z2)));
		
		return - /*1*/ (z2/4)*log((xpp + a)/(-xpp + a)) - /*2*/ (xpp*z2)/(2*a) + /*3*/ (y2/2)*log((xpp + a)/(sqrt(y2 + z2))) + /*4*/ xpp*y2/(2*a) - /*5*/ (ypp*z)*atan((xpp*ypp)/(z*a)) 
					 - /*6*/ ((xpp*y2*z2)/((x2 + z2)*a)) + /*7*/ ((z2*xpp)/(2*a)) + /*8*/ ((xpp*y2*z2)/(2*(x2+z2)*a)) - /*9*/ (xpp*a/2) + /*10*/ (xpp*ypp)*log((ypp + a)/(sqrt(x2+z2))) 
					 - /*11*/ (x3*y2/(2*a*(x2+z2))); 
  break;

  case 1: 

    if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0; 
  	if ((ax<n)&&(az<n)) return - /*6*/(ypp*a/2);
    if ((ay<n)&&(az<n)) return + /*7*/(x2/2)*log((a + ypp)/sqrt(x2+z2)); 
    if (ax<n) return - /*5*/(z2/4)*log((ypp+a)/(-ypp+a)) - /*6*/(ypp*a/2);
  	if (ay<n) return - /*5*/(z2/4)*log((ypp+a)/(-ypp+a)) + /*7*/(x2/2)*log((a + ypp)/sqrt(x2+z2));	
    if (az<n) return + /*2*/(xpp*ypp)*log((xpp + a)/sqrt(y2+z2)) - /*3*/((x2*y3)/(2*a*(y2 + z2))) - /*6*/(ypp*a/2) + /*7*/(x2/2)*log((a + ypp)/sqrt(x2+z2)) + /*8*/(x2*ypp/(2*a));

		return - /*1*/((x2*ypp*z2)/(2*(y2 + z2)*a)) + /*2*/(xpp*ypp)*log((xpp + a)/sqrt(y2+z2)) - /*3*/((x2*y3)/(2*a*(y2 + z2))) - /*4*/(xpp*z)*atan((xpp*ypp)/(z*a))  
					 - /*5*/(z2/4)*log((ypp+a)/(-ypp+a)) - /*6*/(ypp*a/2) + /*7*/(x2/2)*log((a + ypp)/sqrt(x2+z2)) + /*8*/(x2*ypp/(2*a));
	break;

  case 2:

    if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return + /*6*/(z*a/2) + /*7*/z3/(2*a); 
  	if ((ax<n)&&(az<n)) return 0;
    if ((ay<n)&&(az<n)) return 0; 
    if (ax<n) return + /*6*/(z*a/2) + /*7*/z3/(2*a) - /*8*/(ypp*z/2)*log((ypp + a)/(-ypp + a)) + /*9*/(y2*z3)/(2*a*(x2 + z2));
  	if (ay<n) return - /*1*/((xpp*z)/2)*log((xpp + a)/(-xpp + a)) + /*2*/(x2*z3)/(2*a*(y2 + z2)) + /*6*/(z*a/2) + /*7*/z3/(2*a);	
    if (az<n) return 0;

    return - /*1*/((xpp*z)/2)*log((xpp + a)/(-xpp + a)) + /*2*/(x2*z3)/(2*a*(y2 + z2)) - /*3*/(x2*y2*z)/(2*a*(y2 + z2)) - /*4*/(xpp*ypp)*atan(xpp*ypp/(z*a)) 
					 + /*5*/(x2*y2*z)*(x2 + y2 + 2*z2)/((x2 + z2)*(y2 + z2)*a) + /*6*/(z*a/2) + /*7*/z3/(2*a) - /*8*/(ypp*z/2)*log((ypp + a)/(-ypp + a)) + /*9*/(y2*z3)/(2*a*(x2 + z2))  
					 - /*10*/(x2*y2*z/(2*a*(x2 + z2)));
  	break;}
}


/////////////////////////////////////////////////////////
//////////////////  Average volume  /////////////////////
/////////////////////////////////////////////////////////

//scalar potential of volume created by volume, approxiation
double aprox_a(double *r, double *r0)			
{
	double pi=0,rpp[3]={0},mi0=0;

	pi=4.0*atan(1.0);
	mi0=4.0*pi*1e-7;	

	vec_subs2(r,r0,rpp);

return mi0/(4*pi*vec_mod(rpp));		
}

//scalar potential of volume created by volume, q total charge
double spvv(double *r, double *r0, double q, double x[3], double y[3])			
{
	double pi, e,a,b,c,a1,b1,c1;
	double ep=8.854187817620e-12;
	double rpp[3]={0};
	double phi=0, phi1=0;
	double mFVV[4][4][4]={0}, I[3][4]={0};	
	int i, j, k, m, n, o,p;	
	double FVV (double xpp,double ypp, double zpp);
	pi=4.0*atan(1.0);

	a=x[0];	a1=y[0];
	b=x[1]; b1=y[1];
	c=x[2]; c1=y[2];

	vec_subs2(r,r0,rpp);

	e  = (q/(a*b*c*a1*b1*c1))/(4*pi*ep);
	
 	I[0][0] = rpp [0] + (a1/2) + (a/2); //x1
	I[0][1] = rpp [0] + (a1/2) - (a/2);	//x2
	I[0][2] = rpp [0] - (a1/2) + (a/2);	//x3
	I[0][3] = rpp [0] - (a1/2) - (a/2); //x4
	I[1][0] = rpp [1] + (b1/2) + (b/2);	//y1
	I[1][1] = rpp [1] + (b1/2) - (b/2); //y2
	I[1][2] = rpp [1] - (b1/2) + (b/2); //y3
	I[1][3] = rpp [1] - (b1/2) - (b/2); //y4 
	I[2][0] = rpp [2] + (c1/2) + (c/2);	//z1
	I[2][1] = rpp [2] + (c1/2) - (c/2); //z2
	I[2][2] = rpp [2] - (c1/2) + (c/2); //z3
	I[2][3] = rpp [2] - (c1/2) - (c/2); //z4

//loop to evaluate fundamental function FVV with boundary

for(i=0; i<=3; i++)
   { 
			for (j=0; j<=3; j++)
				{ 		
	for (k=0; k<=3; k++)
							{	
									mFVV[i][j][k] = FVV ( I[0][i], I[1][j], I[2][k] );	
									if ((i == 0)||(i == 3)){m=0;}
											else m=1;					
									if ((j == 0)||(j == 3)){n=0;}
											else n=1;
									if ((k == 0)||(k == 3)){o=0;}
											else o=1;

  								mFVV[i][j][k] = pow(-1,(m+n+o))*mFVV[i][j][k];
									phi1+=mFVV[i][j][k];
							}
 				}	
		}
  phi= e * phi1; 
	return phi;
}

//fundamental function scalar potential of volume by volume without boundary
double FVV (double xpp, double ypp, double zpp)
{
	double phi1, x2, x3, x8, y2, y3, z2, x4, x6, y4, y8, y6, z3, z4, z6, a, a3, a5, b, c, d, e, f, n, ax, ay, az;
	x2 = pow(xpp,2.0);  y2 = pow(ypp,2.0);  	z2 = pow(zpp,2.0);
	x3 = pow(xpp,3.0);  y3 = pow(ypp,3.0);    z3 = pow(zpp,3.0);
	x4 = pow(xpp,4.0);  y4 = pow(ypp,4.0);    z4 = pow(zpp,4.0);
	x6 = pow(xpp,6.0);  y6 = pow(ypp,6.0);    z6 = pow(zpp,6.0);
	x8 = pow(xpp,8.0);	y8 = pow(ypp,8.0);

	a  = sqrt(x2 + y2 + z2);
	a3 = pow(a,3.0);
  a5 = pow(a,5.0);
	ax=fabs(xpp);ay=fabs(ypp);az=fabs(zpp);
	n=1e-15; //range of zero point (okolie bodu)

//	cout << "   " << ax << " " << ay << " " << az << " " << endl;
		
  if ((ax<n)&&(ay<n)&&(az<n)) {return 0;}
  if ((ax<n)&&(ay<n)) {//cout << "1 " <<  /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2)) << endl; 
											 return + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2));} 

  if ((ax<n)&&(az<n)) {//cout << "2 " << + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2)) << endl; 
											 return + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2));}

  if ((ay<n)&&(az<n)) {//cout << "3 " << + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2)) << endl; 
											 return + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2));} 

  if (ax<n) {//cout << "4 " << + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2)) + /*8*/(ypp*z4/12)*log(x2 + z2) + /*10*/(y4*zpp/12)*log(x2 + y2)	- /*25*/((ypp*z4)/48)*log(pow(y2 + z2 + ypp*a,2.0) + x2*z2)	- /*26*/((y4*zpp)/48)*log(pow(y2 + z2 + zpp*a,2.0) + x2*y2) <<endl;
										return 	+ /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2)) + /*8*/(ypp*z4/12)*log(x2 + z2) + /*10*/(y4*zpp/12)*log(x2 + y2) 
										- /*25*/((ypp*z4)/48)*log(pow(y2 + z2 + ypp*a,2.0) + x2*z2)	- /*26*/((y4*zpp)/48)*log(pow(y2 + z2 + zpp*a,2.0) + x2*y2);}


  if (ay<n) {//cout <<"5 "<< + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2)) + /*9*/(x4*zpp/8)*log(x2 + y2) + /*12*/(xpp*z4/48)*log(y2 + z2) + /*15*/(x4*zpp/12)*asinh(zpp/sqrt(x2+y2))	- /*23*/((xpp*z4)/48)*log(pow(x2 + z2 + xpp*a,2.0) + y2*z2) - /*24*/((x4*zpp)/16)*log(pow(x2 + z2 + zpp*a,2.0) + x2*y2) <<endl;
						return 	+ /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2)) + /*9*/(x4*zpp/8)*log(x2 + y2) + /*12*/(xpp*z4/48)*log(y2 + z2) + /*15*/(x4*zpp/12)*asinh(zpp/sqrt(x2+y2))
											- /*23*/((xpp*z4)/48)*log(pow(x2 + z2 + xpp*a,2.0) + y2*z2) - /*24*/((x4*zpp)/16)*log(pow(x2 + z2 + zpp*a,2.0) + x2*y2);}



	if (az<n) {//cout << "6 " << endl;cout << - /*22*/((x4*ypp)/16)*log(pow(x2 + y2 + ypp*a,2.0) + x2*z2)<< " " << ((x4*ypp)/16) << " " << pow(x2 + y2 + ypp*a,2.0) << " " << x2 + y2 + (ypp*a) << " " << x2*z2 << 								" " << - /*22*/(x4*ypp)*log(xpp + ypp + ypp*a + xpp*zpp) << " " << xpp + ypp + ypp*a + xpp*zpp << endl; 

					if((x2 + y2 + xpp*a<n)&&(x2 + y2 + ypp*a<n)){
						return  + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2))	+ /*7*/(x4*ypp/8)*log(x2 + z2) + /*11*/(5*xpp*y4/48)*log(y2 + z2) + /*13*/(x4*ypp/12)*asinh(ypp/sqrt(x2+z2)) 
										 + /*14*/(xpp*y4/12)*asinh(xpp/sqrt(y2+z2));}

					if(x2 + y2 + xpp*a<n){//cout << "6a " << + /*14*/(xpp*y4/12)*asinh(xpp/sqrt(y2+z2)) <<  " " << xpp << " " << a <<  " " << (xpp*y4/12)*asinh(xpp/sqrt(y2+z2)) << endl;

						return  + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2))	+ /*7*/(x4*ypp/8)*log(x2 + z2) + /*11*/(5*xpp*y4/48)*log(y2 + z2) + /*13*/(x4*ypp/12)*asinh(ypp/sqrt(x2+z2)) 
										 + /*14*/(xpp*y4/12)*asinh(xpp/sqrt(y2+z2)) - /*22*/((x4*ypp)/16)*log(pow(x2 + y2 + ypp*a,2.0) + x2*z2);}

					if(x2 + y2 + ypp*a<n){
						return  + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2))	+ /*7*/(x4*ypp/8)*log(x2 + z2) + /*11*/(5*xpp*y4/48)*log(y2 + z2) + /*13*/(x4*ypp/12)*asinh(ypp/sqrt(x2+z2)) 
										 + /*14*/(xpp*y4/12)*asinh(xpp/sqrt(y2+z2)) - /*21*/((xpp*y4)/16)*log(pow(x2 + y2 + xpp*a,2.0) + y2*z2);}

						return  + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2))	+ /*7*/(x4*ypp/8)*log(x2 + z2) + /*11*/(5*xpp*y4/48)*log(y2 + z2) + /*13*/(x4*ypp/12)*asinh(ypp/sqrt(x2+z2)) 
										 + /*14*/(xpp*y4/12)*asinh(xpp/sqrt(y2+z2)) - /*21*/((xpp*y4)/16)*log(pow(x2 + y2 + xpp*a,2.0) + y2*z2) - /*22*/((x4*ypp)/16)*log(pow(x2 + y2 + ypp*a,2.0) + x2*z2);} 

///*13*/(x4*ypp/12)*log(ypp + a) == (x4*ypp/12)*asinh(ypp/sqrt(x2+z2))
///*14*/(xpp*y4/12)*log(xpp + a) == (xpp*y4/12)*asinh(xpp/sqrt(y2+z2))
///*15*/(x4*zpp/12)*log(zpp + a) == (x4*zpp/12)*asinh(zpp/sqrt(x2+y2))

	return
				 + /*1*/(a/60)*(pow(x2+y2+z2,2.0) - 5*(x2*y2 + x2*z2 + y2*z2))	

				 + /*2*/(x2*ypp*z2/4)*log(ypp + a) + /*3*/(xpp*y2*z2/4)*log(xpp + a) + /*4*/(x2*y2*zpp/4)*log(zpp + a)

				 - /*5*/(x2*ypp*z2/48) - /*6*/(x2*y2*zpp/48) 

				 + /*7*/(x4*ypp/8)*log(x2 + z2) + /*8*/(ypp*z4/12)*log(x2 + z2) + /*9*/(x4*zpp/8)*log(x2 + y2) + /*10*/(y4*zpp/12)*log(x2 + y2) + /*11*/(5*xpp*y4/48)*log(y2 + z2) + /*12*/(xpp*z4/48)*log(y2 + z2)

				 + /*13*/(x4*ypp/12)*asinh(ypp/sqrt(x2+z2)) + /*14*/(xpp*y4/12)*asinh(xpp/sqrt(y2+z2)) + /*15*/(x4*zpp/12)*asinh(zpp/sqrt(x2 + y2))

				 - /*16*/(xpp*ypp*z3/6)*atan(xpp*ypp/(zpp*a)) - /*17*/(xpp*y3*zpp/6)*atan(xpp*zpp/(ypp*a)) - /*18*/(x3*ypp*zpp/6)*atan(ypp*zpp/(xpp*a)) 
				 + /*19*/(xpp*ypp*z3)*atan(zpp/xpp) + /*20*/(xpp*y3*zpp)*atan(ypp/xpp)

				 - /*21*/((xpp*y4)/16)*log(pow(x2 + y2 + xpp*a,2.0) + y2*z2)
				 - /*22*/((x4*ypp)/16)*log(pow(x2 + y2 + ypp*a,2.0) + x2*z2) 
				 - /*23*/((xpp*z4)/48)*log(pow(x2 + z2 + xpp*a,2.0) + y2*z2) 
				 - /*24*/((x4*zpp)/16)*log(pow(x2 + z2 + zpp*a,2.0) + x2*y2) 
				 - /*25*/((ypp*z4)/48)*log(pow(y2 + z2 + ypp*a,2.0) + x2*z2)
				 - /*26*/((y4*zpp)/48)*log(pow(y2 + z2 + zpp*a,2.0) + x2*y2);

//  if ((ax<n)&&(ay<n)&&(az<n)) {return 0;}
//  if ((ax<n)&&(ay<n)) return (- /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540)); 
//  if ((ax<n)&&(az<n)) return (- /*2*/((y2*a3)/216)- /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540) + /*34*/((5*y4*a)/72));
//  if ((ay<n)&&(az<n)) return (- /*1*/((x2*a3)/216) + /*14*/((5*x4*a)/72) - /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540)); 
//  if (ax<n) return (- /*2*/((y2*a3)/216) + /*10*/((ypp*z4)/16)*log(x2 + z2) + /*12*/((y4*zpp)/16)*log(x2 + y2) 
//										- /*29*/((ypp*z4)/48)*log(((pow((y2 + z2 + (ypp*a)),2.0)) + (x2*z2))/(y4*z6*(x2 + z2))) - /*31*/((y4*zpp)/48)*log(((pow((y2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(y6*z4*(x2 + y2)))
//  	 								- /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540) + /*34*/((5*y4*a)/72) - /*35*/((y2*z2*a)/72));
//  if (ay<n) return (- /*1*/((x2*a3)/216) + /*11*/((x4*zpp)/48)*log(x2 + y2) + /*14*/((5*x4*a)/72) + /*15*/((x4*zpp)/24)*log((zpp + a)/(-zpp + a)) - /*24*/((x2*z2*a)/72)  
//										- /*26*/((xpp*z4)/48)*log(((pow((x2 + z2 + (xpp*a)),2.0)) + (y2*z2))/(x4*z6*(y2 + z2))) - /*27*/((x4*zpp)/12)*log(((pow((x2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(x6*z4*(x2 + y2)))
//         			      + /*30*/((x4*zpp)/48)*log(((pow((x2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(x8*z4*(x2 + y2))) - /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540));
//	if (az<n) return (- /*1*/((x2*a3)/216) - /*2*/((y2*a3)/216) + /*6*/((x2*y2*a)/18) + /*7*/((x4*ypp)/24)*log((ypp + a)/(-ypp + a)) + /*9*/((x4*ypp)/48)*log(x2 + z2)
//										+ /*13*/((xpp*y4)/24)*log((xpp + a)/(-xpp + a)) + /*14*/((5*x4*a)/72) - /*22*/((xpp*y4)/12)*log(((pow((x2 + y2 + (xpp*a)),2.0)) + (y2*z2))/(x4*y6*(y2 + z2)))
//			     		      - /*23*/((x4*ypp)/12)*log(((pow((x2 + y2 + (ypp*a)),2.0)) + (x2*z2))/(x6*y4*(x2 + z2))) + /*25*/((xpp*y4)/48)*log(((pow((x2 + y2 + (xpp*a)),2.0)) + (y2*z2))/(x4*y8*(y2 + z2)))
//			  						+ /*28*/((x4*ypp)/48)*log(((pow((x2 + y2 + (ypp*a)),2.0)) + (x2*z2))/(x8*y4*(x2 + z2))) - /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540) + /*34*/((5*y4*a)/72));

//	return(- /*1*/((x2*a3)/216) - /*2*/((y2*a3)/216) + /*3*/((x2*ypp*z2)/4)*log(ypp + a) + /*4*/((xpp*y2*z2)/4)*log(xpp + a) + /*5*/((x2*y2*zpp)/4)*log(zpp + a)
//			   + /*6*/((x2*y2*a)/18) + /*7*/((x4*ypp)/24)*log((ypp + a)/(-ypp + a)) - /*8*/((x2*ypp*z2)/48) + /*9*/((x4*ypp)/48)*log(x2 + z2) + /*10*/((ypp*z4)/16)*log(x2 + z2) 
//				 + /*11*/((x4*zpp)/48)*log(x2 + y2) + /*12*/((y4*zpp)/16)*log(x2 + y2) + /*13*/((xpp*y4)/24)*log((xpp + a)/(-xpp + a)) + /*14*/((5*x4*a)/72) + /*15*/((x4*zpp)/24)*log((zpp + a)/(-zpp + a))
//				 - /*16*/((x2*y2*zpp)/48) + /*17*/((xpp*ypp*z3)/6)*atan((y2 + z2 + (ypp*a))/(xpp*zpp)) + /*18*/((xpp*y3*zpp)/6)*atan((y2 + z2 + (zpp*a))/(xpp*ypp)) - /*19*/((xpp*ypp*z3)/3)*atan((xpp*ypp)/(zpp*a))
//				 - /*20*/((xpp*y3*zpp)/3)*atan((xpp*zpp)/(ypp*a)) - /*21*/((x3*ypp*zpp)/6)*atan((ypp*zpp)/(xpp*a)) - /*22*/((xpp*y4)/12)*log(((pow((x2 + y2 + (xpp*a)),2.0)) + (y2*z2))/(x4*y6*(y2 + z2)))
//				 - /*23*/((x4*ypp)/12)*log(((pow((x2 + y2 + (ypp*a)),2.0)) + (x2*z2))/(x6*y4*(x2 + z2))) - /*24*/((x2*z2*a)/72) + /*25*/((xpp*y4)/48)*log(((pow((x2 + y2 + (xpp*a)),2.0)) + (y2*z2))/(x4*y8*(y2 + z2)))
//				 - /*26*/((xpp*z4)/48)*log(((pow((x2 + z2 + (xpp*a)),2.0)) + (y2*z2))/(x4*z6*(y2 + z2))) - /*27*/((x4*zpp)/12)*log(((pow((x2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(x6*z4*(x2 + y2))) 
//			   + /*28*/((x4*ypp)/48)*log(((pow((x2 + y2 + (ypp*a)),2.0)) + (x2*z2))/(x8*y4*(x2 + z2))) - /*29*/((ypp*z4)/48)*log(((pow((y2 + z2 + (ypp*a)),2.0)) + (x2*z2))/(y4*z6*(x2 + z2)))
//				 + /*30*/((x4*zpp)/48)*log(((pow((x2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(x8*z4*(x2 + y2))) - /*31*/((y4*zpp)/48)*log(((pow((y2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(y6*z4*(x2 + y2)))
//				 - /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540) + /*34*/((5*y4*a)/72) - /*35*/((y2*z2*a)/72)); 
}

//vector electric field of volume by volume, q total charge 
double vefvv (double *r, double *r0, double q, double a, double b, double c, double a1, double b1, double c1, int l)
{	
	double FEVV(int s, double xpp, double ypp, double zpp);
	double pi,phi=0,d;
	double ep=8.854187817620e-12;
	double rpp[3]={0},I[3][4]={0},mFEVV[4][4][4]={0} ;
	int i, j, k,m, n, o;
	pi=4.0*atan(1.0);
	vec_subs2(r,r0,rpp);

	d  = (q/(a*b*c*a1*b1*c1))/(4*pi*ep);
	I[0][0] = rpp [0] + (a1/2) + (a/2); 	//x1
	I[0][1] = rpp [0] + (a1/2) - (a/2);		//x2
	I[0][2] = rpp [0] - (a1/2) + (a/2);		//x3
	I[0][3] = rpp [0] - (a1/2) - (a/2); 	//x4
	I[1][0] = rpp [1] + (b1/2) + (b/2);		//y1
	I[1][1] = rpp [1] + (b1/2) - (b/2); 	//y2
	I[1][2] = rpp [1] - (b1/2) + (b/2); 	//y3
	I[1][3] = rpp [1] - (b1/2) - (b/2); 	//y4 
	I[2][0] = rpp [2] + (c1/2) + (c/2);		//z1
	I[2][1] = rpp [2] + (c1/2) - (c/2); 	//z2
 	I[2][2] = rpp [2] - (c1/2) + (c/2); 	//z3
	I[2][3] = rpp [2] - (c1/2) - (c/2); 	//z4

//loop to evaluate fundamental function FEVV with boundary

phi=0;
	for(i=0; i<=3; i++)
  	 { 
			for (j=0; j<=3; j++)
		     { 

						for (k=0; k<=3; k++)
							{		mFEVV[i][j][k] = FEVV (l, I[0][i], I[1][j], I[2][k] );	
									if ((i == 0)||(i == 3)){m=0;}
											else m=1;					
									if ((j == 0)||(j == 3)){n=0;}
											else n=1;
									if ((k == 0)||(k == 3)){o=0;}
											else o=1;		
  						  	mFEVV[i][j][k] = pow(-1,(m+n+o))*mFEVV[i][j][k];
								  phi+=mFEVV[i][j][k];
							}
 			    }
	   }
  return d * (-1) * phi; 
}
	
//magnetic field of volume created by volume, approxiation , r0-prime

double Bvv(double *r, double *r0, double V, int i)			
{double pi=0,rpp[3]={0},mi0,p;

	pi=4.0*atan(1.0);
	mi0=4.0*pi*1e-7;	

	vec_subs2(r,r0,rpp);

	p=pow(vec_mod(rpp),3.0);

 	return  (mi0*V*(r[i]-r0[i]))/(4*pi*p);		
}

//magnetic field of volume by volume
double magfvv (double *r, double *r0, double a[3], double b[3], int l)
{	
	double FEVV(int s, double xpp, double ypp, double zpp);
	double pi, mi0, phi=0, d, ep=8.854187817620e-12, rpp[3]={0}, I[3][4]={0}, mFEVV[4][4][4]={0};
	int i, j, k,m, n, o;
	pi=4.0*atan(1.0);
	mi0=4.0*pi*1e-7;
	vec_subs2(r,r0,rpp);

	d  = mi0/(4*pi*b[0]*b[1]*b[2]);

	I[0][0] = rpp [0] + (b[0]/2) + (a[0]/2);	//x1
	I[0][1] = rpp [0] + (b[0]/2) - (a[0]/2);	//x2
	I[0][2] = rpp [0] - (b[0]/2) + (a[0]/2);	//x3
	I[0][3] = rpp [0] - (b[0]/2) - (a[0]/2);	//x4
	I[1][0] = rpp [1] + (b[1]/2) + (a[1]/2);	//y1
	I[1][1] = rpp [1] + (b[1]/2) - (a[1]/2);	//y2
	I[1][2] = rpp [1] - (b[1]/2) + (a[1]/2);	//y3
	I[1][3] = rpp [1] - (b[1]/2) - (a[1]/2);	//y4 
	I[2][0] = rpp [2] + (b[2]/2) + (a[2]/2);	//z1
	I[2][1] = rpp [2] + (b[2]/2) - (a[2]/2);	//z2
 	I[2][2] = rpp [2] - (b[2]/2) + (a[2]/2);	//z3
	I[2][3] = rpp [2] - (b[2]/2) - (a[2]/2);	//z4

//loop to evaluate fundamental function FEVV with boundary

phi=0;
	for(i=0; i<=3; i++)
  	 { 
			for (j=0; j<=3; j++)
		     { 

						for (k=0; k<=3; k++)
							{		mFEVV[i][j][k] = FEVV (l, I[0][i], I[1][j], I[2][k] );	
									if ((i == 0)||(i == 3)){m=0;}
											else m=1;					
									if ((j == 0)||(j == 3)){n=0;}
											else n=1;
									if ((k == 0)||(k == 3)){o=0;}
											else o=1;		
  						  	mFEVV[i][j][k] = pow(-1,(m+n+o))*mFEVV[i][j][k];
								  phi+=mFEVV[i][j][k];
							}
 			    }
	   }
  return d * (-1) * phi; 
}
		
//fundamental function vector electric field of volume by volume without boundary
double FEVV (int l, double xpp, double ypp, double zpp)
{
	double x2,y2,z2,x3,x4,x6,x8,y3,y4,y6,y8,z3,z4,z6,z8,a,b,c,d,phi,ax,ay,az,n;
	x2 = pow(xpp,2.0);  y2 = pow(ypp,2.0);  z2 = pow(zpp,2.0);
	x3 = pow(xpp,3.0);	y3 = pow(ypp,3.0);	z3 = pow(zpp,3.0);
	x4 = pow(xpp,4.0);	y4 = pow(ypp,4.0);	z4 = pow(zpp,4.0);
	x6 = pow(xpp,6.0);	y6 = pow(ypp,6.0);	z6 = pow(zpp,6.0);	
	x8 = pow(xpp,8.0);	y8 = pow(ypp,8.0);	z8 = pow(zpp,8.0);	
	a  = sqrt (x2 + y2 + z2) ;
	ax=fabs(xpp);
	ay=fabs(ypp);
	az=fabs(zpp);
	n=1e-16; //range of zero point (okolie bodu)
switch (l) {
  case 0:
  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return + /*6*/(z4/24)*log((xpp + a)/((-1)*xpp + a)) + /*8*/(z4/48)*log(z2 + y2);
  	if ((ax<n)&&(az<n)) return /*9*/(y4/16)*log(z2 + y2);  
  	if ((ay<n)&&(az<n)) return - /*16*/(xpp/9)*pow(a,3.0) + /*18*/((7*x3*a)/36); 
  	if (ax<n) return /*5*/((z2*y2)/4)*log(xpp + a) + /*6*/(z4/24)*log((xpp + a)/((-1)*xpp + a)) - /*7*/((z2*y2)/48) + /*8*/(z4/48)*log(z2 + y2) + /*9*/(y4/16)*log(z2 + y2) 
										 + /*11*/((zpp*y3)/6)*atan((y2 + x2 + (xpp*a))/(zpp*ypp));
  	if (ay<n) return - /*1*/((z2*xpp*a)/72) + /*6*/(z4/24)*log((xpp + a)/((-1)*xpp + a)) + /*8*/(z4/48)*log(z2 + y2) - /*12*/((zpp*x3)/12)*log((pow((z2 + x2 + (zpp*a)),2.0) + (y2*x2))/(z4*x6*(y2 + x2))) 
				 						 - /*13*/(z4/12)*log((pow((z2 + x2 + (xpp*a)),2.0) + (z2*y2))/(z6*x4*(z2 + y2))) - /*16*/(xpp/9)*pow(a,3.0) + /*18*/((7*x3*a)/36) 
										 + /*21*/(z4/48)*log((pow((z2 + x2 + (xpp*a)),2.0) + (z2*y2))/(z8*x4*(z2 + y2))); 
		if (az<n) return /*3*/((ypp*x3)/4)*log(z2 + x2) + /*9*/(y4/16)*log(z2 + y2) - /*16*/(xpp/9)*pow(a,3.0) + /*18*/((7*x3*a)/36) 
										 - /*19*/((x3*ypp)/12)*log((pow((y2 + x2 + (ypp*a)),2.0) + (z2*x2))/(y4*x6*(z2 + x2))) - /*20*/((y2*xpp*a)/72) 
										 - /*22*/(y4/48)*log((pow((y2 + x2 + (xpp*a)),2.0) + (z2*y2))/(y6*x4*(z2 + y2)));  
	
  	return - /*1*/((z2*xpp*a)/72) + /*2*/((z2*ypp*xpp)/2)*log(ypp + a) + /*3*/((ypp*x3)/4)*log(z2 + x2) + /*4*/((zpp*y2*xpp)/2)*log(zpp + a) + /*5*/((z2*y2)/4)*log(xpp + a) 
 	  	  	 + /*6*/(z4/24)*log((xpp + a)/((-1)*xpp + a)) - /*7*/((z2*y2)/48) + /*8*/(z4/48)*log(z2 + y2) + /*9*/(y4/16)*log(z2 + y2) + /*10*/((zpp*ypp*x2)/2)*atan((y2 + x2 + (ypp*a))/(zpp*xpp)) 
					 + /*11*/((zpp*y3)/6)*atan((y2 + x2 + (xpp*a))/(zpp*ypp)) - /*12*/((zpp*x3)/12)*log((pow((z2 + x2 + (zpp*a)),2.0) + (y2*x2))/(z4*x6*(y2 + x2))) 
					 - /*13*/(z4/12)*log((pow((z2 + x2 + (xpp*a)),2.0) + (z2*y2))/(z6*x4*(z2 + y2))) - /*14*/(zpp*ypp*x2)*atan((zpp*ypp)/(xpp*a)) - /*15*/((2*zpp*y3)/6)*atan((zpp*xpp)/(ypp*a)) 
					 - /*16*/(xpp/9)*pow(a,3.0) - /*17*/((z3*ypp)/6)*atan((ypp*xpp)/(zpp*a)) + /*18*/((7*x3*a)/36) - /*19*/((x3*ypp)/12)*log((pow((y2 + x2 + (ypp*a)),2.0) + (z2*x2))/(y4*x6*(z2 + x2))) 
					 - /*20*/((y2*xpp*a)/72) + /*21*/(z4/48)*log((pow((z2 + x2 + (xpp*a)),2.0) + (z2*y2))/(z8*x4*(z2 + y2))) - /*22*/(y4/48)*log((pow((y2 + x2 + (xpp*a)),2.0) + (z2*y2))/(y6*x4*(z2 + y2))); 
  		break;
	case 1: 
  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return + /*9*/(z4/16)*log(x2 + z2);
  	if ((ax<n)&&(az<n)) return - /*16*/(ypp/9)*pow(a,3.0) + /*18*/((7*y3*a)/36);
  	if ((ay<n)&&(az<n)) return + /*6*/(x4/24)*log((ypp + a)/((-1)*ypp + a)) + /*8*/(x4/48)*log(x2 + z2);
  	if (ax<n) return /*3*/((zpp*y3)/4)*log(x2 + y2) + /*9*/(z4/16)*log(x2 + z2) - /*16*/(ypp/9)*pow(a,3.0) + /*18*/((7*y3*a)/36) 
									 - /*19*/((y3*zpp)/12)*log((pow((z2 + y2 + (zpp*a)),2.0) + (x2*y2))/(z4*y6*(x2 + y2))) - /*20*/((z2*ypp*a)/72) 
		 							 - /*22*/(z4/48)*log((pow((z2 + y2 + (ypp*a)),2.0) + (x2*z2))/(z6*y4*(x2 + z2)));
 	  if (ay<n) return /*5*/((x2*z2)/4)*asinh(ypp/sqrt(x2 + z2)) + /*6*/(x4/24)*log((ypp + a)/((-1)*ypp + a)) - /*7*/((x2*z2)/48) + /*8*/(x4/48)*log(x2 + z2) + /*9*/(z4/16)*log(x2 + z2)
                   + /*11*/((xpp*z3)/6)*atan((z2 + y2 + (ypp*a))/(xpp*zpp));
  	if (az<n) return - /*1*/((x2*ypp*a)/72) + /*6*/(x4/24)*log((ypp + a)/((-1)*ypp + a)) + /*8*/(x4/48)*log(x2 + z2) - /*12*/((xpp*y3)/12)*log((pow((x2 + y2 + (xpp*a)),2.0) + (z2*y2))/(x4*y6*(z2 + y2))) 
										 - /*13*/(x4/12)*log((pow((x2 + y2 + (ypp*a)),2.0) + (x2*z2))/(x6*y4*(x2 + z2))) - /*16*/(ypp/9)*pow(a,3.0) + /*18*/((7*y3*a)/36) 
										 + /*21*/(x4/48)*log((pow((x2 + y2 + (ypp*a)),2.0) + (x2*z2))/(x8*y4*(x2 + z2))); 

		return - /*1*/((x2*ypp*a)/72) + /*2*/((x2*zpp*ypp)/2)*asinh(zpp/sqrt(x2 + y2)) + /*3*/((zpp*y3)/4)*log(x2 + y2) + /*4*/((xpp*z2*ypp)/2)*asinh(xpp/sqrt(z2 + y2)) 
					 + /*5*/((x2*z2)/4)*asinh(ypp/sqrt(x2 + z2)) + /*6*/(x4/24)*log((ypp + a)/((-1)*ypp + a)) - /*7*/((x2*z2)/48) + /*8*/(x4/48)*log(x2 + z2) + /*9*/(z4/16)*log(x2 + z2) 
					 + /*10*/((xpp*zpp*y2)/2)*atan((z2 + y2 + (zpp*a))/(xpp*ypp)) + /*11*/((xpp*z3)/6)*atan((z2 + y2 + (ypp*a))/(xpp*zpp)) 
					 - /*12*/((xpp*y3)/12)*log((pow((x2 + y2 + (xpp*a)),2.0) + (z2*y2))/(x4*y6*(z2 + y2))) - /*13*/(x4/12)*log((pow((x2 + y2 + (ypp*a)),2.0) + (x2*z2))/(x6*y4*(x2 + z2))) 
					 - /*14*/(xpp*zpp*y2)*atan((xpp*zpp)/(ypp*a)) - /*15*/((2*xpp*z3)/6)*atan((xpp*ypp)/(zpp*a)) - /*16*/(ypp/9)*pow(a,3.0) - /*17*/((x3*zpp)/6)*atan((zpp*ypp)/(xpp*a)) + /*18*/((7*y3*a)/36) 
					 - /*19*/((y3*zpp)/12)*log((pow((z2 + y2 + (zpp*a)),2.0) + (x2*y2))/(z4*y6*(x2 + y2))) - /*20*/((z2*ypp*a)/72) + /*21*/(x4/48)*log((pow((x2 + y2 + (ypp*a)),2.0) + (x2*z2))/(x8*y4*(x2 + z2))) 
					 - /*22*/(z4/48)*log((pow((z2 + y2 + (ypp*a)),2.0) + (x2*z2))/(z6*y4*(x2 + z2)));	
			break;
  case 2:
		if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return - /*16*/(zpp/9)*pow(a,3.0) + /*18*/((7*z3*a)/36);
  	if ((ax<n)&&(az<n)) return + /*9*/(y4/16)*log(x2 + y2);  
  	if ((ay<n)&&(az<n)) return /*8*/(x4/48)*log(x2 + y2); 
  	if (ax<n) return /*3*/((ypp*z3)/4)*log(x2 + z2) + /*9*/(y4/16)*log(x2 + y2) - /*16*/(zpp/9)*pow(a,3.0) + /*18*/((7*z3*a)/36) 
										 - /*19*/((z3*ypp)/12)*log((pow((y2 + z2 + (ypp*a)),2.0) + (x2*z2))/(y4*z6*(x2 + z2))) - /*20*/((y2*zpp*a)/72) 
										 - /*22*/(y4/48)*log((pow((y2 + z2 + (zpp*a)),2.0) + (x2*y2))/(y6*z4*(x2 + y2)));
 		if (ay<n) return - /*1*/((x2*zpp*a)/72) + /*6*/(x4/24)*log((zpp + a)/((-1)*zpp + a)) + /*8*/(x4/48)*log(x2 + y2) - /*12*/((xpp*z3)/12)*log((pow((x2 + z2 + (xpp*a)),2.0) + (y2*z2))/(x4*z6*(y2 + z2))) 
    	               - /*13*/(x4/12)*log((pow((x2 + z2 + (zpp*a)),2.0) + (x2*y2))/(x6*z4*(x2 + y2))) - /*16*/(zpp/9)*pow(a,3.0) + /*18*/((7*z3*a)/36) 
										 + /*21*/(x4/48)*log((pow((x2 + z2 + (zpp*a)),2.0) + (x2*y2))/(x8*z4*(x2 + y2)));
		if (az<n) return + /*5*/((x2*y2)/4)*log(zpp + a) - /*7*/((x2*y2)/48) + /*8*/(x4/48)*log(x2 + y2) + /*9*/(y4/16)*log(x2 + y2) + /*11*/((xpp*y3)/6)*atan((y2 + z2 + (zpp*a))/(xpp*ypp)); 

 		return - /*1*/((x2*zpp*a)/72) + /*2*/((x2*ypp*zpp)/2)*log(ypp + a) + /*3*/((ypp*z3)/4)*log(x2 + z2) + /*4*/((xpp*y2*zpp)/2)*log(xpp + a) + /*5*/((x2*y2)/4)*log(zpp + a) 			
					 + /*6*/(x4/24)*log((zpp + a)/((-1)*zpp + a)) - /*7*/((x2*y2)/48) + /*8*/(x4/48)*log(x2 + y2) + /*9*/(y4/16)*log(x2 + y2) + /*10*/((xpp*ypp*z2)/2)*atan((y2 + z2 + (ypp*a))/(xpp*zpp)) 
					 + /*11*/((xpp*y3)/6)*atan((y2 + z2 + (zpp*a))/(xpp*ypp)) - /*12*/((xpp*z3)/12)*log((pow((x2 + z2 + (xpp*a)),2.0) + (y2*z2))/(x4*z6*(y2 + z2))) 
				   - /*13*/(x4/12)*log((pow((x2 + z2 + (zpp*a)),2.0) + (x2*y2))/(x6*z4*(x2 + y2))) - /*14*/(xpp*ypp*z2)*atan((xpp*ypp)/(zpp*a)) - /*15*/((2*xpp*y3)/6)*atan((xpp*zpp)/(ypp*a)) 
					 - /*16*/(zpp/9)*pow(a,3.0) - /*17*/((x3*ypp)/6)*atan((ypp*zpp)/(xpp*a)) + /*18*/((7*z3*a)/36) - /*19*/((z3*ypp)/12)*log((pow((y2 + z2 + (ypp*a)),2.0) + (x2*z2))/(y4*z6*(x2 + z2))) 
					 - /*20*/((y2*zpp*a)/72) + /*21*/(x4/48)*log((pow((x2 + z2 + (zpp*a)),2.0) + (x2*y2))/(x8*z4*(x2 + y2))) - /*22*/(y4/48)*log((pow((y2 + z2 + (zpp*a)),2.0) + (x2*y2))/(y6*z4*(x2 + y2)));  
  	  break;
			}
}

//vector potential of volume created by volume
double vpvv(double *r, double *r0, double a, double b, double c, double a1, double b1, double c1, int l)			
{
	double pi, mi0, e, A1=0, A2=0;
	double rpp[3]={0};
	double mF[4][4][4]={0}, I[3][4]={0};	
	int i, j, k, m, n, o;	
	double FVP4 (int l, double xpp,double ypp, double zpp);

	pi=4.0*atan(1.0);
	mi0=4.0*pi*1e-7;
	vec_subs2(r,r0,rpp);
	e  = (mi0/(4*pi*a*b*c*a1*b1*c1));		

 	I[0][0] = rpp [0] + (a1/2) + (a/2); //x1
	I[0][1] = rpp [0] + (a1/2) - (a/2);	//x2
	I[0][2] = rpp [0] - (a1/2) + (a/2);	//x3
	I[0][3] = rpp [0] - (a1/2) - (a/2); //x4
	I[1][0] = rpp [1] + (b1/2) + (b/2);	//y1
	I[1][1] = rpp [1] + (b1/2) - (b/2); //y2
	I[1][2] = rpp [1] - (b1/2) + (b/2); //y3
	I[1][3] = rpp [1] - (b1/2) - (b/2); //y4 
	I[2][0] = rpp [2] + (c1/2) + (c/2);	//z1
	I[2][1] = rpp [2] + (c1/2) - (c/2); //z2
	I[2][2] = rpp [2] - (c1/2) + (c/2); //z3
	I[2][3] = rpp [2] - (c1/2) - (c/2); //z4

//loop to evaluate fundamental function FVV with boundary

for(i=0; i<=3; i++)
   { 
			for (j=0; j<=3; j++)
				{ 		
						for (k=0; k<=3; k++)
							{	
									mF[i][j][k] = FVP4 ( l, I[0][i], I[1][j], I[2][k] );	
									if ((i == 0)||(i == 3)){m=0;}
											else m=1;					
									if ((j == 0)||(j == 3)){n=0;}
											else n=1;
									if ((k == 0)||(k == 3)){o=0;}
											else o=1;

  								mF[i][j][k] = pow(-1,(m+n+o))*mF[i][j][k];	
									A1+=mF[i][j][k];
							}
 				}	
		}

  return e * (A1 + A2); 
}

//fundamental function vector potential of volume by volume without boundary (old formula scalar potential)
double FVP3 (int l, double xpp, double ypp, double zpp)
{
	double phi1, x2, x3, x8, y2, y3, z2, x4, x6, y4, y8, y6, z3, z4, z6, a, a3, a5, b, c, d, e, f, n, ax, ay, az;
	x2 = pow(xpp,2.0);  y2 = pow(ypp,2.0);  	z2 = pow(zpp,2.0);
	x3 = pow(xpp,3.0);  y3 = pow(ypp,3.0);    z3 = pow(zpp,3.0);
	x4 = pow(xpp,4.0);  y4 = pow(ypp,4.0);    z4 = pow(zpp,4.0);
	x6 = pow(xpp,6.0);  y6 = pow(ypp,6.0);    z6 = pow(zpp,6.0);
	x8 = pow(xpp,8.0);	y8 = pow(ypp,8.0);

	a  = sqrt(x2 + y2 + z2);
	a3 = pow(a,3.0);
  a5 = pow(a,5.0);
	ax=fabs(xpp);ay=fabs(ypp);az=fabs(zpp);
	n=1e-16; //range of zero point (okolie bodu)
switch (l) {
  case 0:
		if ((ax<n)&&(ay<n)&&(az<n)) {return 0;}
		if ((ax<n)&&(ay<n)) return (- /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540)); 
		if ((ax<n)&&(az<n)) return (- /*2*/((y2*a3)/216)- /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540) + /*34*/((5*y4*a)/72));
		if ((ay<n)&&(az<n)) return (- /*1*/((x2*a3)/216) + /*14*/((5*x4*a)/72) - /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540)); 
		if (ax<n) return (- /*2*/((y2*a3)/216) + /*10*/((ypp*z4)/16)*log(x2 + z2) + /*12*/((y4*zpp)/16)*log(x2 + y2) 
											- /*29*/((ypp*z4)/48)*log(((pow((y2 + z2 + (ypp*a)),2.0)) + (x2*z2))/(y4*z6*(x2 + z2))) - /*31*/((y4*zpp)/48)*log(((pow((y2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(y6*z4*(x2 + y2)))
			 								- /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540) + /*34*/((5*y4*a)/72) - /*35*/((y2*z2*a)/72));
		if (ay<n) return (- /*1*/((x2*a3)/216) + /*11*/((x4*zpp)/48)*log(x2 + y2) + /*14*/((5*x4*a)/72) + /*15*/((x4*zpp)/24)*log((zpp + a)/(-zpp + a)) - /*24*/((x2*z2*a)/72)  
											- /*26*/((xpp*z4)/48)*log(((pow((x2 + z2 + (xpp*a)),2.0)) + (y2*z2))/(x4*z6*(y2 + z2))) - /*27*/((x4*zpp)/12)*log(((pow((x2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(x6*z4*(x2 + y2)))
			     			      + /*30*/((x4*zpp)/48)*log(((pow((x2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(x8*z4*(x2 + y2))) - /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540));
		if (az<n) return (- /*1*/((x2*a3)/216) - /*2*/((y2*a3)/216) + /*6*/((x2*y2*a)/18) + /*7*/((x4*ypp)/24)*log((ypp + a)/(-ypp + a)) + /*9*/((x4*ypp)/48)*log(x2 + z2)
											+ /*13*/((xpp*y4)/24)*log((xpp + a)/(-xpp + a)) + /*14*/((5*x4*a)/72) - /*22*/((xpp*y4)/12)*log(((pow((x2 + y2 + (xpp*a)),2.0)) + (y2*z2))/(x4*y6*(y2 + z2)))
						 		      - /*23*/((x4*ypp)/12)*log(((pow((x2 + y2 + (ypp*a)),2.0)) + (x2*z2))/(x6*y4*(x2 + z2))) + /*25*/((xpp*y4)/48)*log(((pow((x2 + y2 + (xpp*a)),2.0)) + (y2*z2))/(x4*y8*(y2 + z2)))
											+ /*28*/((x4*ypp)/48)*log(((pow((x2 + y2 + (ypp*a)),2.0)) + (x2*z2))/(x8*y4*(x2 + z2))) - /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540) + /*34*/((5*y4*a)/72));

		return(- /*1*/((x2*a3)/216) - /*2*/((y2*a3)/216) + /*3*/((x2*ypp*z2)/4)*log(ypp + a) + /*4*/((xpp*y2*z2)/4)*log(xpp + a) + /*5*/((x2*y2*zpp)/4)*log(zpp + a)
					 + /*6*/((x2*y2*a)/18) + /*7*/((x4*ypp)/24)*log((ypp + a)/(-ypp + a)) - /*8*/((x2*ypp*z2)/48) + /*9*/((x4*ypp)/48)*log(x2 + z2) + /*10*/((ypp*z4)/16)*log(x2 + z2) 
					 + /*11*/((x4*zpp)/48)*log(x2 + y2) + /*12*/((y4*zpp)/16)*log(x2 + y2) + /*13*/((xpp*y4)/24)*log((xpp + a)/(-xpp + a)) + /*14*/((5*x4*a)/72) + /*15*/((x4*zpp)/24)*log((zpp + a)/(-zpp + a))
					 - /*16*/((x2*y2*zpp)/48) + /*17*/((xpp*ypp*z3)/6)*atan((y2 + z2 + (ypp*a))/(xpp*zpp)) + /*18*/((xpp*y3*zpp)/6)*atan((y2 + z2 + (zpp*a))/(xpp*ypp)) - /*19*/((xpp*ypp*z3)/3)*atan((xpp*ypp)/(zpp*a))
					 - /*20*/((xpp*y3*zpp)/3)*atan((xpp*zpp)/(ypp*a)) - /*21*/((x3*ypp*zpp)/6)*atan((ypp*zpp)/(xpp*a)) - /*22*/((xpp*y4)/12)*log(((pow((x2 + y2 + (xpp*a)),2.0)) + (y2*z2))/(x4*y6*(y2 + z2)))
					 - /*23*/((x4*ypp)/12)*log(((pow((x2 + y2 + (ypp*a)),2.0)) + (x2*z2))/(x6*y4*(x2 + z2))) - /*24*/((x2*z2*a)/72) + /*25*/((xpp*y4)/48)*log(((pow((x2 + y2 + (xpp*a)),2.0)) + (y2*z2))/(x4*y8*(y2+z2)))
					 - /*26*/((xpp*z4)/48)*log(((pow((x2 + z2 + (xpp*a)),2.0)) + (y2*z2))/(x4*z6*(y2 + z2))) - /*27*/((x4*zpp)/12)*log(((pow((x2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(x6*z4*(x2 + y2))) 
					 + /*28*/((x4*ypp)/48)*log(((pow((x2 + y2 + (ypp*a)),2.0)) + (x2*z2))/(x8*y4*(x2 + z2))) - /*29*/((ypp*z4)/48)*log(((pow((y2 + z2 + (ypp*a)),2.0)) + (x2*z2))/(y4*z6*(x2 + z2)))
	 				 + /*30*/((x4*zpp)/48)*log(((pow((x2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(x8*z4*(x2 + y2))) - /*31*/((y4*zpp)/48)*log(((pow((y2 + z2 + (zpp*a)),2.0)) + (x2*y2))/(y6*z4*(x2 + y2)))
					 - /*32*/(a5/45) + /*33*/((7*a3*( - 2.0*x2 - 2.0*y2 + 3.0*z2))/540) + /*34*/((5*y4*a)/72) - /*35*/((y2*z2*a)/72)); 
		break;

	case 1: 
  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;
  	if ((ay<n)&&(az<n)) return 0;
  	if (ax<n) return 0;
 	  if (ay<n) return 0;
  	if (az<n) return 0; 

		return 0;	
			break;

  case 2:
		if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;  
  	if ((ay<n)&&(az<n)) return 0; 
  	if (ax<n) return 0;
 		if (ay<n) return 0;
		if (az<n) return 0;
 		return 0;  
  	  break;
			}
}

//fundamental function vector potential of volume by volume without boundary
double FVP4 (int l, double xpp, double ypp, double zpp)
{
	double x2,x3,x4,x5,x6,x7,x8,x10,y2,y3,y4,y5,y6,y7,y8,y10,z2,z3,z4,z5,z6,z7,z8,z10,a,a3,b,b1,b2,c,ax,ay,az,n;
	x2 = pow(xpp,2.0);  y2 = pow(ypp,2.0);  z2 = pow(zpp,2.0);
	x3 = pow(xpp,3.0);	y3 = pow(ypp,3.0);	z3 = pow(zpp,3.0);
	x4 = pow(xpp,4.0);	y4 = pow(ypp,4.0);	z4 = pow(zpp,4.0);
	x5 = pow(xpp,5.0);	y5 = pow(ypp,5.0);	z5 = pow(zpp,5.0);
	x6 = pow(xpp,6.0);	y6 = pow(ypp,6.0);	z6 = pow(zpp,6.0);
	x7 = pow(xpp,7.0);	y7 = pow(ypp,7.0);	z7 = pow(zpp,7.0);
	x8 = pow(xpp,8.0);	y8 = pow(ypp,8.0);	z8 = pow(zpp,8.0);
	x10= pow(xpp,10.0);	y10 = pow(ypp,10.0);z10 = pow(zpp,10.0);	
	a  = sqrt (x2 + y2 + z2) ;
	b  = pow((x2 + y2 + xpp*a),2.0);
	b1  = pow((x2 + z2 + xpp*a),2.0);
	b2  = pow((x2 + z2 + zpp*a),2.0);
  a3 = pow(a,3.0);
	ax=fabs(xpp);
	ay=fabs(ypp);
	az=fabs(zpp);
	n=1e-16; //range of zero point (okolie bodu)
switch (l) {
  case 0:
  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;  
  	if ((ay<n)&&(az<n)) return 0; 
  	if (ax<n) return 0;
  	if (ay<n) return 0; 
		if (az<n) return 0;  
	
  	return - /*1*/(xpp*pow(a,5.0)/60) - /*2*/(17*x3*a3/600) + /*3*/(47*xpp*y2*a3/900) - /*4*/(101*xpp*z2*a3/1800) - /*5*/(17*xpp*y4*a/400) + /*6*/(11*y6/480)*log(xpp + a) + /*7*/(y4*z2/16)*log(xpp + a)
					 + /*8*/(83*x5*a/5400) - /*9*/(409*x3*y2*a/5400) + /*10*/(13*x3*z2*a/2700) + /*11*/(37*xpp*y2*z2*a/1800) - /*12*/(13*xpp*z4*a/3600) + /*13*/(3*y6/160)*log(-xpp + a) - /*14*/(z6/72)*log(xpp + a)
					 + /*15*/(y2*z4/48)*log(xpp + a) + /*16*/(xpp*y4*zpp/24)*log(zpp + a)  - /*17*/(x4*ypp*zpp/12)*atan(ypp*zpp/(xpp*a)) + /*18*/(x5*ypp/24)*log(-ypp + a) + /*19*/(x5*zpp/24)*log(-zpp + a) 
					 - /*20*/(37*x3*ypp*z2/240) - /*21*/(xpp*ypp*z4/160) - /*22*/(ypp*z5/60)*atan((y2 + z2 + ypp*a)/(xpp*zpp)) + /*23*/(120*y2*log(400*(b + y2*z2)/(x4*y10*(y2 +z2)))) 
					 + /*24*/(120*x5*log(400*(b + x2*z2)/(x10*y4*(x2 +z2)))) + /*25*/(z6/144)*log(-xpp + a) - /*26*/(y5*zpp/60)*atan((y2 + z2 + zpp*a)/(xpp*ypp)) + /*27*/(y6/720)*log(36*(b + y2*z2)/(x2*y6*(y2 +z2)))
					 + /*28*/(z6/720)*log(400*(b1 + y2*z2)/(x4*z10*(y2 +z2)))	+ /*29*/(x5*zpp/120)*log(400*(b2 + x2*y2)/(x10*z4*(x2 +y2))) + /*30*/(x3*ypp*z2/12)*log(ypp + a) + /*31*/(xpp*ypp*z4/24)*log(ypp + a)
					 + /*32*/(x3*y2*zpp/12)*log(zpp + a);
  		break;
	case 1: 
  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;
  	if ((ay<n)&&(az<n)) return 0;
  	if (ax<n) return 0;
 	  if (ay<n) return 0;
  	if (az<n) return 0; 
		return 0;	
			break;
  case 2:
		if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;  
  	if ((ay<n)&&(az<n)) return 0; 
  	if (ax<n) return 0;
 		if (ay<n) return 0;
		if (az<n) return 0;
 		return 0;  
  	  break;
			}
}

//fundamental function vector potential of volume by volume without boundary
double FVP5 (int l, double xpp, double ypp, double zpp)
{
	double x2,x3,x4,x5,x6,x7,x8,x10,x12,y2,y3,y4,y5,y6,y7,y8,y10,z2,z3,z4,z5,z6,z7,z8,z10,a,a3,a5,b,b1,b2,b3,b4,b5,c,ax,ay,az,n;
	x2 = pow(xpp,2.0);  y2 = pow(ypp,2.0);  z2 = pow(zpp,2.0);
	x3 = pow(xpp,3.0);	y3 = pow(ypp,3.0);	z3 = pow(zpp,3.0);
	x4 = pow(xpp,4.0);	y4 = pow(ypp,4.0);	z4 = pow(zpp,4.0);
	x5 = pow(xpp,5.0);	y5 = pow(ypp,5.0);	z5 = pow(zpp,5.0);
	x6 = pow(xpp,6.0);	y6 = pow(ypp,6.0);	z6 = pow(zpp,6.0);
	x7 = pow(xpp,7.0);	y7 = pow(ypp,7.0);	z7 = pow(zpp,7.0);
	x8 = pow(xpp,8.0);	y8 = pow(ypp,8.0);	z8 = pow(zpp,8.0);
	x10= pow(xpp,10.0);	y10 = pow(ypp,10.0);z10 = pow(zpp,10.0);	
	x12= pow(xpp,12.0);
	a  = sqrt (x2 + y2 + z2) ;
	b  = pow((x2 + y2 + xpp*a),2.0);
	b1  = pow((x2 + z2 + xpp*a),2.0);
	b2  = pow((x2 + z2 + zpp*a),2.0);
	b3  = pow((x2 + y2 + ypp*a),2.0);
	b4  = pow((y2 + z2 + ypp*a),2.0);
	b5  = pow((y2 + z2 + zpp*a),2.0);
  a3 = pow(a,3.0);
	a5 = pow(a,5.0);
	ax=fabs(xpp);
	ay=fabs(ypp);
	az=fabs(zpp);
	n=1e-16; //range of zero point (okolie bodu)
switch (l) {
  case 0:
  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;  
  	if ((ay<n)&&(az<n)) return 0; 
  	if (ax<n) return 0;
  	if (ay<n) return 0; 
		if (az<n) return 0;  
	
  	return /*1*/(pow(a,7.0)/315) + /*2*/(29*x2*a5/6300) - /*3*/(19*y2*a5/12600) - /*4*/(95*z2*a5/63000) + /*5*/(46*x2*y2*a3/3375) - /*6*/(937*x4*a3/54000) - /*7*/(9*x2*z2*a3/500) 
				 + /*8*/(((41/23625)+(1/48))*x6*a) + /*9*/(7*x6*zpp/1440)*log(zpp + a) + /*10*/(((322752/189000)-(11/400))*x2*y4*a) + /*11*/(x2*y4*zpp/48)*log(zpp + a) + /*12*/(23*y4*a3/54000) 
				 - /*13*/(3*y2*z2*a3/13500) - /*14*/(((41/23625)+(7/3600))*y6*a) + /*15*/(((322752/189000)-(23/900))*x4*y2*a) + /*16*/(164*x4*z2*a/189000) - /*17*/(23*y4*z2*a/94500) 
				 + /*18*/(((328/189000)+(1/144))*x2*y2*z2*a) - /*19*/(123*x2*z4*a/189000) + /*20*/(187*y2*z4*a/126000) + /*21*/(615*z6*a/189000) + /*22*/(53*x6*zpp/1440)*log(-zpp + a) 
				 - /*23*/(y6*zpp/144)*log(zpp + a) - /*24*/(x5*ypp*zpp/15)*atan(ypp*zpp/(xpp*a)) + /*25*/(7*x6*ypp/1440)*log(ypp + a) + /*26*/(53*x6*ypp/1440)*log(-ypp + a)  
				 + /*27*/(x6*ypp/180)*log(36*(b3 + x2*z2)/(x12*y4*(x2 + z2))) + /*28*/(ypp*z6/180)*log(400*(b4 + x2*z2)/(y4*z10*(x2 + z2))) + /*29*/(x6*zpp/180)*log(36*(b2 + x2*y2)/(x12*z4*(x2 + y2)))
				 + /*30*/(y6*zpp/180)*log(400*(b5 + x2*y2)/(y10*z4*(x2 + y2))) - /*31*/(5*x4*ypp*z2/288) - /*32*/(x2*ypp*z4/144) - /*33*/(ypp*z6/144)*log(ypp + a) + /*34*/(x2*ypp*z4/48)*log(ypp + a)
				 + /*35*/(x4*y2*zpp/16)*log(zpp + a);    
					 
  		break;
	case 1: 
  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;
  	if ((ay<n)&&(az<n)) return 0;
  	if (ax<n) return 0;
 	  if (ay<n) return 0;
  	if (az<n) return 0; 

		return 0;	
			break;
  case 2:
		if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;  
  	if ((ay<n)&&(az<n)) return 0; 
  	if (ax<n) return 0;
 		if (ay<n) return 0;
		if (az<n) return 0;
 		return 0;  
  	  break;
			}
}

//fundamental function vector potential of volume by volume without boundary
double FVP6 (int l, double xpp, double ypp, double zpp)
{
	double x2,x3,x4,x5,x6,x7,x8,x10,x12,y2,y3,y4,y5,y6,y7,y8,y10,y12,z2,z3,z4,z5,z6,z7,z8,z10,a,a3,a5,b,b1,b2,b3,b4,b5,b6,c,ax,ay,az,n;
	x2 = pow(xpp,2.0);  y2 = pow(ypp,2.0);  z2 = pow(zpp,2.0);
	x3 = pow(xpp,3.0);	y3 = pow(ypp,3.0);	z3 = pow(zpp,3.0);
	x4 = pow(xpp,4.0);	y4 = pow(ypp,4.0);	z4 = pow(zpp,4.0);
	x5 = pow(xpp,5.0);	y5 = pow(ypp,5.0);	z5 = pow(zpp,5.0);
	x6 = pow(xpp,6.0);	y6 = pow(ypp,6.0);	z6 = pow(zpp,6.0);
	x7 = pow(xpp,7.0);	y7 = pow(ypp,7.0);	z7 = pow(zpp,7.0);
	x8 = pow(xpp,8.0);	y8 = pow(ypp,8.0);	z8 = pow(zpp,8.0);
	x10= pow(xpp,10.0);	y10= pow(ypp,10.0);z10 = pow(zpp,10.0);	
	x12= pow(xpp,12.0);	y12= pow(ypp,12.0);
	a  = sqrt (x2 + y2 + z2) ;
	b  = pow((x2 + y2 + xpp*a),2.0);
	b1  = pow((x2 + z2 + xpp*a),2.0);
	b2  = pow((x2 + z2 + zpp*a),2.0);
	b3  = pow((x2 + y2 + ypp*a),2.0);
	b4  = pow((y2 + z2 + ypp*a),2.0);
	b5  = pow((y2 + z2 + zpp*a),2.0);
	b6  = pow((x2 + y2 + xpp*a),2.0);
  a3 = pow(a,3.0);
	a5 = pow(a,5.0);
	ax=fabs(xpp);
	ay=fabs(ypp);
	az=fabs(zpp);
	n=1e-16; //range of zero point (okolie bodu)
switch (l) {
  case 0:
  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;  
  	if ((ay<n)&&(az<n)) return 0; 
  	if (ax<n) return 0;
  	if (ay<n) return 0; 
		if (az<n) return 0;  
	
  	return - /*1*/(x2*a5/450) + /*2*/(y2*a5/900) + /*3*/(z2*a5/900) - /*4*/(5*x4*ypp*z2/72) + /*5*/(((1/675)-(984/283500))*x4*y2*a) + /*6*/(13*x6*ypp/240)*log(ypp + a) + /*7*/(x4*ypp*z2/8)*log(ypp + a)
					 - /*8*/(881*x4*a3/40500)	+ /*9*/(((1/27)-(328/283500))*x6*a) + /*10*/(((-7/5400)-(328/283500))*y6*a) - /*11*/(2603*x2*y4*a/94500) + /*12*/(((164/283500)-(1/216))*x4*z2*a) 
					 + /*13*/(((328/283500)-(29/1200))*x2*y2*z2*a) + /*14*/(((328/283500)-(1/1350))*y4*z2*a) + /*15*/((123/283500)*x2*z4*a) + /*16*/(((123/283500)+(1/1800))*y2*z4*a) + /*17*/(615*z6*a/283500)
					 - /*18*/(x2*ypp*z4/96) + /*19*/(17*x6*ypp/240)*log(-ypp + a) - /*20*/(ypp*z6/24)*log(ypp + a) - /*21*/(x3*y2*z2/16) + /*22*/(5*x3*y4/144)*log(xpp + a) + /*23*/(x3*z4/48)*log(xpp + a) 
					 + /*24*/(13*x6*zpp/240)*log(zpp + a) + /*25*/(x4*y2*zpp/8)*log(zpp + a) + /*26*/(23*y4*a3/81000) - /*27*/(89*x2*y2*a3/81000) - /*28*/(y2*z2*a3/6750) + /*29*/(13*x2*z2*a3/13500)
					 + /*30*/(17*x6*zpp/240)*log(-zpp + a) - /*31*/(y6*zpp/24)*log(zpp + a) + /*32*/(x3*ypp*z3/18)*atan((y2 + z2 + ypp*a)/(xpp*zpp)) - /*33*/(x3*y4/36)*log(144*(b6 + y2*z2)/(x4*y6*(y2 + z2)))
					 - /*34*/(x6*ypp/36)*log(144*(b3 + x2*z2)/(x6*y4*(x2 + z2))) + /*35*/(x3*y4/144)*log(-xpp + a) + /*36*/(x3*z4/48)*log(-xpp + a) + /*37*/(x6*ypp/45)*log(36*(b3 + x2*z2)/(x12*y4*(x2 + z2)))
					 + /*38*/(ypp*z6/72)*log(144*(b4 + x2*z2)/(y4*z6*(x2 + z2))) + /*39*/(x3*y3*zpp/18)*atan((y2 + z2 + zpp*a)/(xpp*ypp)) + /*40*/(x3*y4/144)*log(16*(b6 + y2*z2)/(x4*y8*(y2 + z2)))
					 - /*41*/(x3*z4/144)*log(144*(b1 + y2*z2)/(x4*z6*(y2 + z2))) - /*42*/(x6*zpp/36)*log(144*(b2 + x2*y2)/(x6*z4*(x2 + y2))) + /*43*/(x6*zpp/45)*log(36*(b2 + x2*y2)/(x12*z4*(x2 + y2)))
					 + /*44*/(y6*zpp/72)*log(144*(b5 + x2*y2)/(y6*z4*(x2 + y2))) - /*45*/(x3*ypp*z3/9)*atan((xpp*ypp)/(zpp*a)) - /*46*/(x3*y3*zpp/9)*atan((xpp*zpp)/(ypp*a)) 
					 - /*47*/(x5*ypp*zpp/10)*atan((ypp*zpp)/(xpp*a)) + /*48*/(26*x2*y4*a/900) + /*49*/(ypp*z6/120)*log(400*(b4 + x2*z2)/(y4*z10*(x2 + z2))) ;
					 + /*50*/(y6*zpp/120)*log(400*(b5 + x2*y2)/(y10*z4*(x2 + y2)));					 

  		break;
	case 1: 
  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;
  	if ((ay<n)&&(az<n)) return 0;
  	if (ax<n) return 0;
 	  if (ay<n) return 0;
  	if (az<n) return 0; 

		return 0;	
			break;
  case 2:
		if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;  
  	if ((ay<n)&&(az<n)) return 0; 
  	if (ax<n) return 0;
 		if (ay<n) return 0;
		if (az<n) return 0;
 		return 0;  
  	  break;
			}
}


//fundamental function vector potential of volume by volume without boundary
double FVP7 (int l, double xpp, double ypp, double zpp)
{
	double x2,x3,x4,x5,x6,x7,x8,x10,x12,y2,y3,y4,y5,y6,y7,y8,y10,y12,z2,z3,z4,z5,z6,z7,z8,z10,a,a3,a5,b,b1,b2,b3,b4,b5,b6,c,ax,ay,az,n;
	x2 = pow(xpp,2.0);  y2 = pow(ypp,2.0);  z2 = pow(zpp,2.0);
	x3 = pow(xpp,3.0);	y3 = pow(ypp,3.0);	z3 = pow(zpp,3.0);
	x4 = pow(xpp,4.0);	y4 = pow(ypp,4.0);	z4 = pow(zpp,4.0);
	x5 = pow(xpp,5.0);	y5 = pow(ypp,5.0);	z5 = pow(zpp,5.0);
	x6 = pow(xpp,6.0);	y6 = pow(ypp,6.0);	z6 = pow(zpp,6.0);
	x7 = pow(xpp,7.0);	y7 = pow(ypp,7.0);	z7 = pow(zpp,7.0);
	x8 = pow(xpp,8.0);	y8 = pow(ypp,8.0);	z8 = pow(zpp,8.0);
	x10= pow(xpp,10.0);	y10= pow(ypp,10.0);z10 = pow(zpp,10.0);	
	x12= pow(xpp,12.0);	y12= pow(ypp,12.0);
	a  = sqrt (x2 + y2 + z2) ;
	b  = pow((x2 + y2 + xpp*a),2.0);
	b1  = pow((x2 + z2 + xpp*a),2.0);
	b2  = pow((x2 + z2 + zpp*a),2.0);
	b3  = pow((x2 + y2 + ypp*a),2.0);
	b4  = pow((y2 + z2 + ypp*a),2.0);
	b5  = pow((y2 + z2 + zpp*a),2.0);
	b6  = pow((x2 + y2 + xpp*a),2.0);
  a3 = pow(a,3.0);
	a5 = pow(a,5.0);
	ax=fabs(xpp);
	ay=fabs(ypp);
	az=fabs(zpp);
	n=1e-16; //range of zero point (okolie bodu)
switch (l) {
  case 0:
  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;  
  	if ((ay<n)&&(az<n)) return 0; 
  	if (ax<n) return 0;
  	if (ay<n) return 0; 
		if (az<n) return 0;  
	
  	return - /*1*/(xpp*a5/360) - /*2*/(593*x3*y2*a/3600) + /*3*/(19*xpp*y4*a/2400) - /*4*/(xpp*y2*z2*a/1200) + /*5*/(ypp*z5/15)*atan(xpp*ypp/(zpp*a)) + /*6*/(y6/240)*log(400*(b6 + y2*z2)/(x4*y10*(z2 + y2)))
					 - /*7*/(x5*ypp/30)*log(400*(b3 + x2*z2)/(x10*y4*(x2 + z2))) - /*8*/(23*x3*z2*a/3600) - /*9*/(xpp*z4*a/2400) - /*10*/(17*y6/2880)*log(xpp + a) + /*11*/(17*y6/2880)*log(-xpp + a) 
					 + /*12*/(z6/288)*log(-xpp + a) + /*13*/(1721*x3*a3/5400) - /*14*/(11*xpp*y2*a3/3600) + /*15*/(19*xpp*z2*a3/3600) + /*16*/(y4*z2/96)*log(xpp + a) + /*17*/(y2*z4/96)*log(xpp + a)
					 + /*18*/(x2*y4/8)*log((xpp + a)/(sqrt(y2 + z2))) + /*19*/(x2*y2*z2/8)*log((xpp + a)/(sqrt(y2 + z2))) - /*20*/(7*x2*y4/96)*log(xpp + a) + /*21*/(x2*z4/32)*log(xpp + a) + /*22*/(x2*y4/32)*log(y2+z2)
				   - /*23*/(x2*z4/32)*log(y2 + z2) + /*24*/(y5*zpp/15)*atan(xpp*zpp/(ypp*a)) + /*25*/(y6/1440)*log(36*(b + y2*z2)/(x4*y12*(y2 + z2))) - /*26*/(z6/180)*log(400*(b1 + y2*z2)/(y4*z5*(y2 + z2)))
					 + /*27*/(x5*zpp/240)*log(400*(b2 + x2*y2)/(x10*z4*(x2 + y2))) + /*28*/(x5*zpp/48)*log(zpp + a) + /*29*/(x3*y2*zpp/6)*log((zpp + a)/(sqrt(x2 + y2))) - /*30*/(19*x3*ypp*z2/240)
 					 + /*31*/(x2*ypp*z3/12)*atan((y2 + z2 + ypp*a)/xpp*zpp) - /*32*/(x2*y4/24)*log(144*(b6 + y2*z2)/(y4*y6*(y2 + z2))) - /*33*/(x5*ypp/24)*log(144*(b3 + x2*z2)/(x6*y6*(x2 + z2)))
					 + /*34*/(x2*y4/96)*log(-xpp + a) + /*35*/(x2*y4/96)*log(16*(b + y2*z2)/(x4*y8*(y2 + z2))) - /*36*/(x2*z4/96)*log(144*(b1 + y2*z2)/(x4*z6*(y2 + z2))) 
					 - /*37*/(x5*zpp/24)*log(144*(b2 + x2*y2)/(x6*z4*(x2 + y2))) - /*38*/(3*y5*zpp/40)*atan((y2 + z2 + zpp*a)/xpp*ypp) + /*39*/(z6/160)*log(400*(b1 + y2*z2)/(x4*z10*(y2 + z2))) 
					 - /*40*/(x2*ypp*z3/6)*atan(xpp*ypp/(zpp*a)) - /*41*/(xpp*ypp*z4/40) + /*42*/(ypp*z5/10)*atan(xpp*zpp/(y2 + z2 + ypp*a)) - /*43*/(x2*y3*zpp/6)*atan(xpp*zpp/(ypp*a)) 
					 - /*44*/(x4*ypp*zpp/8)*atan(ypp*zpp/(xpp*a)) - /*45*/(x5*ypp/16)*log(ypp + a) - /*46*/(x5*ypp/16)*log(-ypp + a) + /*47*/(x5*zpp/16)*log(-zpp + a);
 					 
  		break;
	case 1: 
  	if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;
  	if ((ay<n)&&(az<n)) return 0;
  	if (ax<n) return 0;
 	  if (ay<n) return 0;
  	if (az<n) return 0; 

		return 0;	
			break;
  case 2:
		if ((ax<n)&&(ay<n)&&(az<n)) return 0;
  	if ((ax<n)&&(ay<n)) return 0;
  	if ((ax<n)&&(az<n)) return 0;  
  	if ((ay<n)&&(az<n)) return 0; 
  	if (ax<n) return 0;
 		if (ay<n) return 0;
		if (az<n) return 0;
 		return 0;  
  	  break;
			}
}




