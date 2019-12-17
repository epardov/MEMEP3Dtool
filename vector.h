using namespace std;

//contains several common vector operations for 3-dimensional vectors

///////////////////////////////////////////////////////////////////////////////////////////////////

//Returns the sign of a scalar (double precision).

int sign (double x)
{
	return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//sets a vector to zero

void vec_zero(double *v)
{
	int i;
	for (i=0;i<3;i++) v[i]=0.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//sets absolute value of vector

void vec_fabs(double *u)
{
	int i;
	for (i=0;i<3;i++) u[i]=fabs(u[i]);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Assigns one vector to the other (u:=v).

void vec_cpy(double *u,double *v)
{
	u[0]=v[0];
	u[1]=v[1];
	u[2]=v[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Assigns one vector to the other (u:=v).

void vec_cpy_1(double *u,double *v)
{
	u[0]=v[3];
	u[1]=v[4];
	u[2]=v[5];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Assigns one vector to the other (u:=v).

void vec_cpy_2(double *u,double *v, double d)
{
	u[0]=v[0];
	u[1]=d;
	u[2]=v[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Assigns one vector[n] to the other (u:=v).

void vec_cpy_n(int n, double *u,double *v)
{ int i;
	for (i=0; i<n; i++)	
	u[i]=v[i];

}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Assigns one vector to the other (u:=v).

void vec_cpy_i(int *u,int *v)
{
	u[0]=v[0];
	u[1]=v[1];
	u[2]=v[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//adition one vector to the other (u+=v).

void addition(double *u,double *v)
{
	u[0]+=v[0];
	u[1]+=v[1];
	u[2]+=v[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Adds two vectors and returns the value in u.

void vec_add(double *u,double *v)
{
	int i;
	for (i=0;i<3;i++) u[i]=u[i]+v[i];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Adds two vectors and returns the value in u.

void vec_add2(double *u,double *v, double *w)
{
	int i;
	for (i=0;i<3;i++) w[i]=u[i]+v[i];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Adds two vectors and returns the value in u.

void vec_add3(double *u, double v)
{
	int i;
	for (i=0;i<3;i++) u[i]=u[i]+v;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Substracts two vectors and returns the value in u, u <- u-v.

void vec_subs(double *u,double *v)
{
	int i;
	for (i=0;i<3;i++) u[i]=u[i]-v[i];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Substracts v from u and returns to w: w=u-v.

void vec_subs2(double *u,double *v,double *w)
{
	int i;
	for (i=0;i<3;i++) w[i]=u[i]-v[i];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Modulus of a vector

double vec_mod(double *u)
{
	return sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Modulus square of a vector

double vec_mod2(double *u)
{
	return u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Calculates the vector distance

double vec_mod3(double *u,double *v)
{int i;
 double w[3]={0};

	for (i=0;i<3;i++) {w[i]=u[i]-v[i];}
	return sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Calculates the vector distance

double vec_mod4(double *u,double *v)
{int i;
 double w[3]={0};

	for (i=0;i<3;i++) {w[i]=u[i]-v[i];}
	return pow(sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]),3);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Calculates the vector distance

double vec_mod5(double *u,double *v)
{int i;
 double w[3]={0};

	for (i=0;i<3;i++) {w[i]=u[i]-v[i];}
	return pow(sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]),5);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Scale a vector and return to u

void vec_scale(double *u,double lam)
{
	u[0]=lam*u[0];
	u[1]=lam*u[1];
	u[2]=lam*u[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Scale a vector and return to w

void vec_scale1(double *u,double lam,double *w)
{
	w[0]=lam*u[0];
	w[1]=lam*u[1];
	w[2]=lam*u[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Division a vector and return to w

void vec_divi(double *u,double lam, double *w)
{
	w[0]=u[0]/lam;
	w[1]=u[1]/lam;
	w[2]=u[2]/lam;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Division a vector and return to w

void vec_div(double *w, double *u, int *v)
{
	w[0]=u[0]/v[0];
	w[1]=u[1]/v[1];
	w[2]=u[2]/v[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Scale a vector and return to u

void vec_scale_n(int n, double *u,double lam)
{	int i;
	for (i=0; i<n; i++)
	u[i]=lam*u[i];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Calculates the dot product 

double vec_dot(double *u,double *v)
{
	return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Calculates the vector product(rotor) and return to w, w=u \times v.

void vec_rot(double *u,double *v,double *w)
{
	w[0]=u[1]*v[2]-u[2]*v[1];
	w[1]=u[2]*v[0]-u[0]*v[2];
	w[2]=u[0]*v[1]-u[1]*v[0];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Calculates the vector product(rotor) and return to w, than modulus |W|

double vec_rot1(double *u,double *v)
{ double w[3]={0};
	w[0]=u[1]*v[2]-u[2]*v[1];
	w[1]=u[2]*v[0]-u[0]*v[2];
	w[2]=u[0]*v[1]-u[1]*v[0];

	return sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//Calculates the vector product(rotor) and return to wx,wy,wz

void vec_rot2(double *u,double *v, double& wx, double& wy, double& wz)
{
	wx = u[1]*v[2] - u[2]*v[1];
	wy = u[2]*v[0] - u[0]*v[2];
	wz = u[0]*v[1] - u[1]*v[0];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//calculates cubic volume from sizes of cell saved in vector 

double vol(double *u)
{
	return u[0]*u[1]*u[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//calculates u + i*v

double oper1(double *w,double *u, int *i,double *v)
{
	w[0] = u[0] + i[0]*v[0];
	w[1] = u[1] + i[1]*v[1];
	w[2] = u[2] + i[2]*v[2];	
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//calculates u - v/2  + d/2 

double oper2(double *w,double *u, double *v, double *d, int i)
{
	w[0] = u[0] - v[0]/2.0 + d[0]/2.0; if(i==0)w[0] = u[0] - v[0];
	w[1] = u[1] - v[1]/2.0 + d[1]/2.0; if(i==1)w[1] = u[1] - v[1];
	w[2] = u[2] - v[2]/2.0 + d[2]/2.0; if(i==2)w[2] = u[2] - v[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//calculates u - v/2

double oper3(double *w, double *u, double *v, int s)
{
	w[0] = u[0] - v[0]/2.0; if(s==0) {w[0] = u[0] - v[0];}
	w[1] = u[1] - v[1]/2.0;	if(s==1) {w[1] = u[1] - v[1];}
	w[2] = u[2] - v[2]/2.0;	if(s==2) {w[2] = u[2] - v[2];}
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//calculates u + i*d + d/2

double oper4(double *w, double *u, int *i, double *d, double a,int s)
{
	w[0] = u[0] + i[0]*d[0] + d[0]/2.0; if(s==0) {w[0] = u[0] + i[0]*d[0];}
	w[1] = u[1] + i[1]*a 		+ a/2.0;		if(s==1) {w[1] = u[1] + i[1]*a;		}
	w[2] = u[2] + i[2]*d[2] + d[2]/2.0;	if(s==2) {w[2] = u[2] + i[2]*d[2];}
}

///////////////////////////////////////////////////////////////////////////////////////////////////

//integral sqrt(r^2+(dr/dfi)^2) dfi

double l_spiral_element(double r1, double b, double fi1, double fi2)
{

return  (r1+b*fi2)*sqrt(b*b+pow(r1+b*fi2,2.0))/(2*b) + 0.5*b*log(r1+b*fi2+sqrt(b*b+pow(r1+b*fi2,2.0)))
		 -( (r1+b*fi1)*sqrt(b*b+pow(r1+b*fi1,2.0))/(2*b) + 0.5*b*log(r1+b*fi1+sqrt(b*b+pow(r1+b*fi1,2.0))));
}

///////////////////////////////////////////////////////////////////////////////////////////////////

