/*Ray tracing - A. D. Pierce, section 8.1 (Equation 8-1.10a - 8-1.10b)*/ 
/*The system of ODE is solved using Runge-Kutta Method of order 4*/

//---------------Standard C librabries are included here-------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "mex.h"
//-----------------------Some Static Values--------------------------------------------

#define yposition_for_ray 0.5
#define PI 3.141592

double *vxd, *vyd;
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
//void mexFunction( int nlhs, mxArray *plhs[], 
//		  int nrhs, const mxArray *prhs[] )
int main ()
{

	void Ray_Acoustics(double *xmin, double *xmax, double *ymin, 
            double *ymax,double *tmaxd,int itp,int yn, int tn,int Typ);
	double *xmin, *xmax, *ymin, *ymax,*tmaxd,*Itype,*Yn;
	int itp,Ny,Nt,Typ;
  	//mwSize mrows,ncols;
  
  /* Assign pointers to each input and output. */
  /*	xmin = mxGetPr(prhs[0]);
  	xmax = mxGetPr(prhs[1]);
  	ymin = mxGetPr(prhs[2]);
  	ymax = mxGetPr(prhs[3]);
    vxd  = mxGetPr(prhs[4]);
    vyd  = mxGetPr(prhs[5]);
    tmaxd = mxGetPr(prhs[6]);
    Itype = mxGetPr(prhs[7]);
    itp = *Itype;
    Yn = mxGetPr(prhs[8]);
    Ny= *Yn;
    Nt= *Yn;
    Itype = mxGetPr(prhs[9]);
    Typ=*Itype;*/
    
    if((xmin=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    if((xmax=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    if((ymin=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    if((ymax=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    if((tmaxd=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    
    if((vxd=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    if((vyd=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    
    *xmin=0;
    *xmax=1.0;
    *ymin=-0.5;
    *ymax=1.5;
    *vxd=0;
    *vyd=0.6;
    *tmaxd=0.35;
    itp=1;
    Ny=40000;
    Nt=40000;
    Typ=0;
    
    Ray_Acoustics(xmin,xmax,ymin,ymax,tmaxd,itp,Ny,Nt,Typ);
//	return;    
}
//-------------------------Main program starts here------------------------------------
void Ray_Acoustics (double *xmini, double *xmaxi, double *ymini, 
        double *ymaxi,double *tmaxd,int Itp,int yn, int tn, int Typ) {

//---------------------------Variables are declared here-------------------------------
	int check;
	double x,y,yv,tv,s1[2],s2[2],k1,k2,k3,k4,x1,y1,yray;
	double delY,delt;
    double xmin,xmax,ymin,ymax,tmax;

//-------------------------Functions Declaration---------------------------------------
	double c (double x, double y,int I);
	double vx (double x, double y,int Typ);
	
	double rhs_s1(double x, double y,double s1,double s2,int I,int Typ);
	double rhs_s2(double x, double y,double s1,double s2,int I,int Typ);
	double rhs_x(double x, double y,double s1,double s2,int I,int Typ);
	double rhs_y(double x, double y,double s1,double s2,int I,int Typ);

//---------------------------File pointer to store the output--------------------------
	FILE *fpw,*fpr;
//-------------------------------------------------------------------------------------
//-------------------------Discretizing [0,1]x[0,1]------------------------------------
    xmin=*xmini;xmax=*xmaxi;ymin=*ymini;ymax=*ymaxi;
	delY=(ymax-ymin)/(yn-1); 
    tmax=*tmaxd;delt=tmax/tn;
//-------------------------------------------------------------------------------------
//--------------------------Opening a file for storing the output----------------------
	fpw=fopen("Data/Linear_Wavefront.dat","w");
//-------------------------------------------------------------------------------------
	yv=ymin;check=0;

	while(yv<=ymax)
	{
//---------------------Initialization--------------------------------------------------
		tv = delt;
		y=yv;x=0.0;
		s1[0]=1.0/(c(x,y,Itp)+vx(x,y,Typ));
		s2[0]=0.0;
//-------------------------------------------------------------------------------------
		//printf("Starting of the ray is at (x,y) = (0,%lf)\n",yv);
		if(yv-delY<yposition_for_ray && yv>=yposition_for_ray)
		{
			fpr=fopen("Data/Linear_Ray.dat","w");
			check=1;
			yray=yv;
		}
		while (tv<tmax)
		{
			k1=delt*rhs_s1(x,y,s1[0],s2[0],Itp,Typ);
			k2=delt*rhs_s1(x,y,s1[0]+0.5*k1,s2[0]+0.5*k1,Itp,Typ);
			k3=delt*rhs_s1(x,y,s1[0]+0.5*k2,s2[0]+0.5*k2,Itp,Typ);
			k4=delt*rhs_s1(x,y,s1[0]+0.5*k3,s2[0]+0.5*k3,Itp,Typ);
			s1[1] = s1[0]+0.166666667*(k1+2*(k2+k3)+k4);
	
			k1=delt*rhs_s2(x,y,s1[0],s2[0],Itp,Typ);
			k2=delt*rhs_s2(x,y,s1[0]+0.5*k1,s2[0]+0.5*k1,Itp,Typ);
			k3=delt*rhs_s2(x,y,s1[0]+0.5*k2,s2[0]+0.5*k2,Itp,Typ);
			k4=delt*rhs_s2(x,y,s1[0]+0.5*k3,s2[0]+0.5*k3,Itp,Typ);
			s2[1] = s2[0]+0.166666667*(k1+2*(k2+k3)+k4);

			k1=delt*rhs_x(x,y,s1[0],s2[0],Itp,Typ);
			k2=delt*rhs_x(x+0.5*k1,y+0.5*k1,s1[0],s2[0],Itp,Typ);
			k3=delt*rhs_x(x+0.5*k2,y+0.5*k2,s1[0],s2[0],Itp,Typ);
			k4=delt*rhs_x(x+k3,y+k3,s1[0],s2[0],Itp,Typ);
			x1 = x+0.166666667*(k1+2*(k2+k3)+k4);

			k1=delt*rhs_y(x,y,s1[0],s2[0],Itp,Typ);
			k2=delt*rhs_y(x+0.5*k1,y+0.5*k1,s1[0],s2[0],Itp,Typ);
			k3=delt*rhs_y(x+0.5*k2,y+0.5*k2,s1[0],s2[0],Itp,Typ);
			k4=delt*rhs_y(x+k3,y+k3,s1[0],s2[0],Itp,Typ);
			y1 = y+0.166666667*(k1+2*(k2+k3)+k4);

			x=x1;y=y1;
			s1[0]=s1[1];s2[0]=s2[1];
			tv+=delt;
			if(check==1)
				fprintf(fpr,"%lf %lf\n",x,y);
		}

		if(check==1)
		{
			check=0;
			fclose(fpr);
		}

		fprintf(fpw,"%lf %lf\n",x,y);

		yv+=delY;
	}
//-------------------------------------------------------------------------------------
//--------------------------Closing the file for storing the output----------------------
	fclose(fpw);
//-------------------------------------------------------------------------------------
	printf("The program is completed successfully.\n");
	printf("The wavefront at time = %lf is stored in the file named Linear_Wavefront.dat\n",tmax);
	printf("The ray at y = %lf is stored in the file named Linear_Ray.dat\n",yray);


}

double c (double x, double y,int I)
{
	double speed;
  
    //concave speed function
    if(I==1)
    {
        if (y<0.0) 
            y=0.0;
        else if (y>1.0)
            y=1.0;
        return(2.0-0.5*(cos(PI*(y-0.5))*cos(PI*(y-0.5))));
    }
    else if(I==2)
    {
    //speed with cavity
        speed = 2.5-0.4*y;
        speed = speed - 0.8*exp(-10.0*(y-0.3)*(y-0.3))*exp(-20.0*(x-0.25)*(x-0.25));
        return(speed);		
    }
    return(0);
}

double c_x (double x, double y,int I)
{
    double speed;
	
   if(I==2)
   {
        speed = -0.8*exp(-10.0*(y-0.3)*(y-0.3))*exp(-20.0*(x-0.25)*(x-0.25));
        speed = -40.0*speed*(x-0.25);
   }
   else if(I==1)
       speed = 0;
   return(speed);
}
double c_y (double x, double y,int I)
{
	double speed;
	
    if(I==2)
    {
        speed = -0.8*exp(-10.0*(y-0.3)*(y-0.3))*exp(-20.0*(x-0.25)*(x-0.25));
        speed = -0.4-20.0*speed*(y-0.3);
        return(speed);
    }
    else if(I==1)
    {
        if (y<0 || y>1)
            return(0.0);
        else
            return(PI*cos(PI*(y-0.5))*sin(PI*(y-0.5)));
    }
    return(0);
}
double vy (double x, double y,int T)
{
    if (T==0)
        return(*vyd);
    else
        return(0.5*cos(T*PI*y));
    
}
double vx (double x, double y,int T)
{
    if (T==0)
        return(*vxd);
    else
        return(-0.5*sin(T*PI*x));        
	
/*    return((x>0.05)?0.25*cos(PI*y):0);
	if(y<0.5)
		return(0.5);
	else
		return(-0.5);
 */
}
double vy_x (double x, double y,int T)
{
	return(0.0);
}
double vy_y (double x, double y,int T)
{
    return(-0.5*T*PI*sin(T*PI*y));
}
double vx_x (double x, double y,int T)
{
    
    return(-0.5*T*PI*cos(T*PI*x));
    
}
double vx_y (double x, double y,int T)
{
	return(0);
//	return((x>0.05)?-0.25*PI*sin(PI*y):0);
}
double rhs_s1 (double x, double y,double s1,double s2,int I,int Typ)
{
	double dummy;
	double c(double x,double y,int I);
	double c_x(double x, double y,int I);
	double vx(double x, double y,int Typ);
	double vy(double x, double y,int Typ);
	double vx_x(double x, double y,int Typ);
	double vy_x(double x, double y,int Typ);

	dummy = -1*(1.0-vx(x,y,Typ)*s1-vy(x,y,Typ)*s2)/c(x,y,I);
	//return(dummy*c_x(x,y,I));
	return(dummy*c_x(x,y,I)-s1*vx_x(x,y,Typ)-s2*vy_x(x,y,Typ));
}
double rhs_s2 (double x, double y,double s1,double s2,int I,int Typ)
{
	double dummy;
	double c(double x,double y,int I);
	double c_y(double x, double y,int I);
	double vx(double x, double y,int Typ);
	double vy(double x, double y,int Typ);
	double vx_x(double x, double y,int Typ);
	double vy_x(double x, double y,int Typ);

	dummy = -1*(1.0-vx(x,y,Typ)*s1-vy(x,y,Typ)*s2)/c(x,y,I);
	return(dummy*c_y(x,y,I)-s1*vx_y(x,y,Typ)-s2*vy_y(x,y,Typ));
}
double rhs_x (double x, double y,double s1,double s2,int I,int Typ)
{
	double dummy;
	double c(double x,double y,int I);
	double vx(double x, double y,int Typ);
	double vy(double x, double y,int Typ);

	dummy = (1.0-vx(x,y,Typ)*s1-vy(x,y,Typ)*s2);
	return((c(x,y,I)*c(x,y,I)*s1/dummy)+vx(x,y,Typ));
}
double rhs_y (double x, double y,double s1,double s2,int I,int Typ)
{
	double dummy;
	double c(double x,double y,int I);
	double vx(double x, double y,int Typ);
	double vy(double x, double y,int Typ);

	dummy = (1.0-vx(x,y,Typ)*s1-vy(x,y,Typ)*s2);
	return((c(x,y,I)*c(x,y,I)*s2/dummy)+vy(x,y,Typ));
}
