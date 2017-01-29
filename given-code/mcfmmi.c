/**
 * This file is mainly to create the input files for the initial conditions and the information about the nodes and the elements.
 * Those input files then further act as an input to the main solver functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592

/* ITs = 1  ==>  concave sound speed
 * ITs = 2  ==>  convex sound speed
 * ITs = 3  ==>  cavity sound speed
 */

/* ITw = 1  =>>  normal wavefront from left to right
 * ITw = -1 =>>  normal wavefront from right to left
 */

float F(float x, float y,int Ind)
{
	float f,Y,X;

    if(Ind==1)
    {
        if (y<0.0) 
        	Y=0.0;
        else if (y>1.0)
            Y=1.0;
        else
        	Y=y;

    	f=2.0-0.5*(cos(PI*(Y-0.5))*cos(PI*(Y-0.5)));    
    }
    else if(Ind==2)
    {
        if (y<0.0) 
        	Y=0.0;
        else if (y>1.0)
            Y=1.0;
        else
        	Y=y;

    	f=2.0-0.5*(sin(PI*(Y-0.5))*sin(PI*(Y-0.5)));    
    }
    else if(Ind==10)
    {
        if (x<0.0) 
        	X=0.0;
        else if (x>1.0)
            X=1.0;
        else
        	X=x;

    	f=2.0-0.5*(cos(PI*(X-0.5))*cos(PI*(X-0.5)));    
    }
    else if(Ind==3)
    {
        f = 2.5 - 0.4*y;
        f = f - 0.8*exp(-10.0*(y-0.3)*(y-0.3))*exp(-20.0*(x-0.25)*(x-0.25));
	}
    else if(Ind==4)
    {
        f = 2.5 - 0.4*y;
        f = f - 0.8*exp(-10.0*(y-0.9)*(y-0.9))*exp(-20.0*(x-0.55)*(x-0.55));
	}
    else
        f = 1.0;
    
    return(f);
}

/**
 * This is the function to generate the information about the initial wave front.
 * For the time being, let us consider a planar wavefront which is parallel to the left boundary
 * The different ITW's list down the different types of initial wavefronts that have been used.
 */
void IWF(FILE *fp,int xn,int Yn, double xmin, double ymin,
        float Dx, float Dy,float vxd,float vyd,int ITs,int ITv,int ITw,int N[2])
{
    FILE *fp1;
    
    
    float T,x,y;
    int j,i,Ng,I,J;
    float F(float x, float y,int Ind);
	float vx(float v,float x,float y,int Ind);
	float vy(float v,float x,float y,int Ind);
    double dt(double F0,double delX,double delY,double xmin,
            double xmax,double ymin,double ymax,int i,int j);
    
    N[0]=0; N[1]=0;
    if(ITw==1)                      //normal from left to right
    {
        for(j=1;j<=Yn;j++)
        {
        	fprintf(fp,"1 %d 0\n",j);
            N[0]++;
        }
       
        for(j=1;j<=Yn;j++)
        {
            i=2;
        	x=Dx;
        	y=ymin+(j-1)*Dy;
        	T=F(x,y,ITs)+vx(vxd,x,y,ITv);
        	T=Dx/T;
        	fprintf(fp,"%d %d %f\n",i,j,T);
            N[1]++;
        }
    }
    if(ITw==11)//normal wavefront from left to right & from right to left
    {
        for(j=1;j<=Yn;j++)
        {
            fprintf(fp,"1 %d 0\n",j);
            N[0]++;
        } 
        for(j=1;j<=Yn;j++)
        {
            fprintf(fp,"%d %d 0\n",xn,j);
            N[0]++;
        } 
        for(j=1;j<=Yn;j++)
        {
            i=2;
            x=Dx;
        	y=ymin+(j-1)*Dy;
        	T=F(x,y,ITs)+vx(vxd,x,y,ITv);
        	T=Dx/T;
        	fprintf(fp,"%d %d %f\n",i,j,T);
            N[1]++;
        }
        for(j=1;j<=Yn;j++)
        {
            i=xn-1;
        	x=Dx;
        	y=ymin+(j-1)*Dy;
        	T=F(x,y,ITs)-vx(vxd,x,y,ITv);
        	T=Dx/T;
        	fprintf(fp,"%d %d %f\n",i,j,T);
            N[1]++;
        }
    }
    else if(ITw==10)                      //normal from bottom to top
    {
        for(i=1;i<=xn;i++)
        {
        	fprintf(fp,"%d 1 0\n",i);
            N[0]++;
        }
    
        for(i=1;i<=xn;i++)
        {
            j=2;
        	y=Dy;
        	x=xmin+(i-1)*Dx;
        	T=F(x,y,ITs)+vy(vyd,x,y,ITv);
        	T=Dy/T;
        	fprintf(fp,"%d %d %f\n",i,j,T);
            N[1]++;
        }
    }
    if(ITw==-1)                      //normal from right to left
    {
        for(j=1;j<=Yn;j++)
        {
        	fprintf(fp,"%d %d 0\n",xn,j);
            N[0]++;
        }
    
        for(j=1;j<=Yn;j++)
        {
            i=xn-1;
        	x=Dx;
        	y=ymin+(j-1)*Dy;
        	T=F(x,y,ITs)-vx(vxd,x,y,ITv);
        	T=Dx/T;
        	fprintf(fp,"%d %d %f\n",i,j,T);
            N[1]++;
        }
    }
    else if(ITw==-10)                      //normal from top to bottom
    {
        for(i=1;i<=xn;i++)
        {
        	fprintf(fp,"%d %d 0\n",i,Yn);
            N[0]++;
        }
    
        for(i=1;i<=xn;i++)
        {
            j=Yn-1;
        	y=Dy;
        	x=xmin+(i-1)*Dx;
        	T=F(x,y,ITs)-vy(vyd,x,y,ITv);
        	T=Dy/T;
        	fprintf(fp,"%d %d %f\n",i,j,T);
            N[1]++;
        }
    }
    else if(ITw==5) //Point source
    {
        //Alive points initialization for Point source
        Ng=3;
        
        for(j=1;j<Ng;j++)
        {
            for(i=1;i<Ng;i++)
            {
                fprintf(fp,"%d %d 0\n",i,j);
                N[0]++;
            }
        }
        
        //Narrow Band initialization for Point source
       
        i=Ng;j=2;
        x=Dx;
    	y=ymin+(j-1)*Dy;
    	T=F(x,y,ITs)+vx(vxd,x,y,ITv);
    	T=Dx/T;
    	fprintf(fp,"%d %d %f\n",i,j,T);
        N[1]++;

        j=Ng;i=2;
        y=Dy;
        x=xmin+(i-1)*Dx;
        T=F(x,y,ITs)+vx(vxd,x,y,ITv);
        T=Dy/T;
        fprintf(fp,"%d %d %f\n",i,j,T);
        N[1]++;
        
        j=Ng;i=Ng;
        x=Dx;y=Dy;
        y=ymin+(j-1)*Dy;
        x=xmin+(i-1)*Dx;
        T=F(x,y,ITs)+vx(vxd,x,y,ITv);
        T=sqrt(Dx*Dx+Dy*Dy)/T;
        fprintf(fp,"%d %d %f\n",i,j,T);
        N[1]++;

    }
    else if(ITw==55) //Point source of both sides with interaction
    {
        //Alive points initialization for point source
        Ng=3;
        for(j=1;j<Ng;j++)
        {
            for(i=1;i<Ng;i++)
            {
                fprintf(fp,"%d %d 0\n",i,j);
                N[0]++;
            }
        }
        for(j=Yn-Ng+2;j<Yn+1;j++)
        {
            for(i=xn-Ng+2;i<xn+1;i++)
            {
                fprintf(fp,"%d %d 0\n",i,j);
                N[0]++;
            }
        }
        
        //Narrow Band initialization for source
        
        for(j=1;j<Ng;j++)
        {
            i=Ng;
            x=Dx;
        	y=ymin+(j-1)*Dy;
        	T=F(x,y,ITs)+vx(vxd,x,y,ITv);
        	T=Dx/T;
        	fprintf(fp,"%d %d %f\n",i,j,T);
            N[1]++;
        }
        
        for(i=1;i<Ng;i++)
        {
        	j=Ng;
            y=Dy;
        	x=xmin+(i-1)*Dx;
        	T=F(x,y,ITs)+vy(vyd,x,y,ITv);
        	T=Dy/T;
        	fprintf(fp,"%d %d %f\n",i,j,T);
            N[1]++;
        }
        
        
        
        //Narrow Band initialization for source
        
        for(j=Yn-Ng+2;j<Yn+1;j++)
        {
            i=xn-Ng+1;
            x=Dx;
        	y=ymin+(j-1)*Dy;
        	T=F(x,y,ITs)-vx(vxd,x,y,ITv);
        	T=Dx/T;
        	fprintf(fp,"%d %d %f\n",i,j,T);
            N[1]++;
        }
        
        for(i=xn-Ng+2;i<xn+1;i++)
        {
        	j=Yn-Ng+1;
            y=Dy;
        	x=xmin+(i-1)*Dx;
        	T=F(x,y,ITs)-vy(vyd,x,y,ITv);
        	T=Dy/T;
        	fprintf(fp,"%d %d %f\n",i,j,T);
            N[1]++;
        }
         
    }
    
    else if(ITw==555) // point source at (0,0) and (1,1) in domain [0,1]x[0,1], taking analytic solution at four alive points at (0,0) and four alive points at (1,1) and NB points are taken as the immediate neighbors of the four alive points at (0.0) and (1,1) and are assigned the analytic solution
    {

        //Alive Points (0,0) and (1,1)
        fprintf(fp,"%d %d 0\n",i,j);
        fprintf(fp,"%d %d\n",xn-1,Yn-1);
        N[0]=2;

	
    }
    else if(ITw==99) // Oblique wavefront
    {
//Alive Points
        
        for(j=Yn;j>0;j--)
        {
        	i=1;
            x=i*Dx;y=j*Dy;
            T=F(x,y,ITs);
            T=dt(T,Dx,Dy,xmin,2,ymin,2,i,j);
           // fprintf(fp,"%d %d %f\n",i,j,T);
        	fprintf(fp,"%d %d %f\n",i,j,(Yn-j)*Dy*cos(PI/4.0)/F(x,y,ITs));
            printf("%d %d %f\n",i,j,T);
            N[0]++;
        }
        for(i=2;i<xn+1;i++)
        {
        	j=Yn;
            x=i*Dx;y=j*Dy;
            T=F(x,y,ITs);
            T=dt(T,Dx,Dy,xmin,2,ymin,2,i,j);
            //fprintf(fp,"%d %d %f\n",i,j,T);
        	fprintf(fp,"%d %d %f\n",i,j,i*Dx*cos(PI/4.0)/F(x,y,ITs));
            N[0]++;
        }
//Narrow Band initialization for source
        fp1=fopen("NBN.dat","r");
        fscanf(fp1,"%d",&I);
        N[1]=I;
        fclose(fp1);
        
        fp1=fopen("NB.dat","r");
        for(j=0;j<N[1];j++)
        {
                fscanf(fp1,"%d%d%f",&I,&J,&T);
                fprintf(fp,"%d %d %f\n",I,J,T);
        }
        fclose(fp1);
//         for(j=Yn-1;j>0;j--)
//         {
//         	i=2;
//             x=i*Dx;y=j*Dy;
//         	fprintf(fp,"%d %d %f\n",i,j,(Yn-j)*Dy*cos(PI/4.0)/F(x,y,ITs));
//             N[1]++;
//         }
//         for(i=3;i<xn+1;i++)
//         {
//         	j=Yn-1;
//             x=i*Dx;y=j*Dy;
//         	fprintf(fp,"%d %d %f\n",i,j,i*Dx*cos(PI/4.0)/F(x,y,ITs));
//             N[1]++;
//         }
	
    }
    return;
}

double dt(double F0,double delX,double delY,double xmin,double xmax,
        double ymin,double ymax,int i,int j)
{
	double distance,c,x,y,m;	
    double modl(double z0);
    
	x=xmin+(i-1)*delX;
	y=ymin+(j-1)*delY;
	m=1.0; // slope of line
	c=ymax-xmin; // found by assuming that the line passes through (xmin,ymax)

	distance=(modl(m*x-y+c))/(sqrt(m*m+1));
	return(distance/F0);
}
double modl(double z0)
{
	if(z0>=0)
		return(z0);
	else
		return(-1*z0);
}
float AnalyticPS(float F,float vx,float vy,float x,float y,float x0,float y0,float T)
{
    float norm,velX,sp;
    
    velX = vx*(x-x0) + vy*(y-y0);
    norm = ((x-x0)*(x-x0)) + ((y-y0)*(y-y0));
    sp   = (F*F - (vx*vx+vy*vy));
    
    T = T + (norm/(velX + sqrt(sp*norm + velX*velX)));
    return(T);
}
float vx(float v,float x,float y,int Ind)
{
	if(Ind==0)
        return(v);
    else
        return(0.5*sin(Ind*PI*x));
}
float vy(float v,float x,float y,int Ind)
{
    if(Ind==0)
        return(v);
    else
        return(0.25*cos(Ind*PI*y));
}


//void mexFunction( int nlhs, mxArray *plhs[], 
//		  int nrhs, const mxArray *prhs[] )
int main ()
{

	double *xmin, *xmax, *ymin, *ymax;
    double *xnp,*ynp,*wind,*sind,*vind;
    double *vy,*vx;
	int xnpt,ynpt,indw,inds,indv;
    void finput (double xmin, double xmax, double ymin, double ymax,
                double vx, double vy, int xn, int Yn, int ITs, int ITv,int ITw);
  /*
	xmin = mxGetPr(prhs[0]);
  	xmax = mxGetPr(prhs[1]);
  	ymin = mxGetPr(prhs[2]);
  	ymax = mxGetPr(prhs[3]);
    
    vx = mxGetPr(prhs[4]);
    vy = mxGetPr(prhs[5]);
  	
    xnp = mxGetPr(prhs[6]);
	xnpt=*xnp;
	ynp = mxGetPr(prhs[7]);
	ynpt=*ynp;
    sind = mxGetPr(prhs[8]);
	inds=*sind;
    vind = mxGetPr(prhs[9]);
	indv=*vind;
    wind = mxGetPr(prhs[10]);
	indw=*wind;*/
    
    if((xmin=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    if((xmax=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    if((ymin=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    if((ymax=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    
    if((vx=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    if((vy=(double *) malloc (3*sizeof(double)))==NULL)
        return(0);
    
    *xmin=0;
    *xmax=1.1;
    *ymin=-0.5;
    *ymax=1.5;
    *vx=0;
    *vy=0.6;

    xnpt=1600;
    ynpt=xnpt;
    inds=1;
    indv=0;
    indw=1;
           
	
    finput(*xmin,*xmax,*ymin,*ymax,*vx,*vy,xnpt,ynpt,inds,indv,indw);
    
	return(0);    
}

void finput (double xmin, double xmax, double ymin, double ymax,
                double vxd, double vyd, int xn, int Yn, int ITs, int ITv, int ITw) 
{
 
    FILE *fp;
    int i,j,NAP,NNP,N[2];
    int I[xn],J[Yn],NI[xn],NJ[Yn];
    float T;
    float Dx, Dy;
    double x,y;
    float F(float x, float y,int Ind);
	float vx(float v,float x,float y,int Ind);
	float vy(float v,float x,float y,int Ind);
    void IWF(FILE *fp,int xn,int Yn, double xmin, double ymin,
        float Dx, float Dy,float vxd,float vyd,int ITs,int ITv,int ITw,int N[2]);

    
    Dx=(xmax-xmin)/(xn-1);
    Dy=(ymax-ymin)/(Yn-1);
    
    fp=fopen("Data/inputD.dat","w");
    
    y=ymin;
    for(j=0;j<Yn;j++,y+=Dy)
    {
    	x=xmin;
     	for(i=0;i<xn;i++,x+=Dx)
         	fprintf(fp,"%f %f %f\n",F(x,y,ITs),vx(vxd,x,y,ITv),vy(vyd,x,y,ITv));
    }
    IWF(fp,xn,Yn,xmin,ymin,Dx,Dy,vxd,vyd,ITs,ITv,ITw,N); 
    fclose(fp);
    
    NAP=N[0];NNP=N[1];
    fp=fopen("Data/inputN.dat","w");
    fprintf(fp,"%d %d %d %d\n",xn,Yn,NAP,NNP);
    fprintf(fp,"%f %f\n",Dx,Dy);
    fclose(fp);
   
    printf("Initial condition is generated successfully\n");
    
    fp=fopen("Data/F.dat","w");
    y=ymin;
	for(j=0;j<Yn;j++,y+=Dy)
    {
    	x=xmin;
     	for(i=0;i<xn;i++,x+=Dx)
         	fprintf(fp,"%f ",F(x,y,ITs));
        fprintf(fp,"\n");
    }
	fclose(fp);
}
