#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Far_pt	 -10000
#define Alive_pt -20000
#define Infinity 1000000.0
#define PI 3.141592

struct sv/// Stores all the information about the solution
{
	double **T;/// stores the minimum wave approach time
	int **I;/// stores the state of the node, whether it is Far away, Alive or Narrow-band
	int **Q;
	int **NIt;
	int **NB;
    int NBP[1];
};

struct mv// This struct contains the information about the given parameters
{
	float **F,**vx,**vy;
	/// F is the given wave speed at all the points in the grid
	/// vx is the speed of the medium in the x-direction
	/// vy is the speed of the medium in the y-direction
};

//-------------------------------------------------------------------------------------
int mmalloc(struct mv *m,int NGx,int NGy)/// Assigning the space to uder-defined parameters
{
	int i;

    NGx=NGx+2;
    NGy=NGy+2;

	if((m->F=(float **) malloc (NGy*sizeof(float)))==NULL)/// Assigning wave speed variable space to store the adresses
        return(0);
    for(i=0;i<NGy;i++)/// Assigning each of the addresses the size of a float
    {
       if((*(m->F+i)=(float *) malloc (NGx*sizeof(float)))==NULL)
           return(0);
    }
	/// Same routine for vx and vy
	if((m->vx=(float **) malloc (NGy*sizeof(float)))==NULL)
        return(0);
    for(i=0;i<NGy;i++)
    {
        if((*(m->vx+i)=(float *) malloc (NGx*sizeof(float)))==NULL)
            return(0);
    }
	if((m->vy=(float **) malloc (NGy*sizeof(float)))==NULL)
        return(0);
    for(i=0;i<NGy;i++)
    {
    	if((*(m->vy+i)=(float *) malloc (NGx*sizeof(float)))==NULL)
            return(0);
    }
	/// Returning 1 if everyting goes fine, for checking purposes
	return(1);
}
//-------------------------------------------------------------------------------------
int smalloc(struct sv *s,int NGx,int NGy,int NNP)/// Assigning the memory to the soltuions parameters
{
	int i,j;

    NGx=NGx+2;
    NGy=NGy+2;

	/// Fpr all the parameter assigning the same way as done for mmalloc
	if((s->T=(double **) malloc (NGy*sizeof(double)))==NULL)
        return(0);
    for(i=0;i<NGy;i++)
    {
        if((*(s->T+i)=(double *) malloc (NGx*sizeof(double)))==NULL)
            return(0);
    }
	if((s->I=(int **) malloc (2*NGy*sizeof(int)))==NULL)
    	return(0);
    for(i=0;i<NGy;i++)
    {
    	if((*(s->I+i)=(int *) malloc (2*NGx*sizeof(int)))==NULL)
            return(0);
    }
    if((s->Q=(int **) malloc (2*NGy*sizeof(int)))==NULL)
    	return(0);
    for(i=0;i<NGy;i++)
    {
    	if((*(s->Q+i)=(int *) malloc (2*NGx*sizeof(int)))==NULL)
            return(0);
    }
    if((s->NIt=(int **) malloc (2*NGy*sizeof(int)))==NULL)
    	return(0);
    for(i=0;i<NGy;i++)
    {
    	if((*(s->NIt+i)=(int *) malloc (2*NGx*sizeof(int)))==NULL)
            return(0);
    }
	/// Intially labelling all points as Far away points and consequently setting the wave approack time to be infinite
//--------------------------------Initialization --------------------------------------
    for(j=0;j<NGy;j++)
	{
		for(i=0;i<NGx;i++)
		{
			*(*(s->I+j)+i)=Far_pt;
			*(*(s->T+j)+i)=Infinity;
		}
	}
//-------------------------------------------------------------------------------------

    if((s->NB=(int **) malloc (2*sizeof(int)))==NULL)
    	return(0);
    for(i=0;i<2;i++)
    {
    	if((*(s->NB+i)=(int *) malloc (2*NNP*sizeof(int)))==NULL)
        	return(0);
    }
	return(1);
}
//------------------------------------------------------------------------------------
/// Function that finds(and returns) the minimum time of all solutions in the narrow band
double minim (int **NB,int No_of_NB_Pts,int **Ind_Type,double  **T)
{
	int i,j,k,dummy;
    double Td;

    i=No_of_NB_Pts-1;
    Td=T[NB[0][i]][NB[1][i]];
    for(j=No_of_NB_Pts-2;j>=0;j--)
    {
        if(T[NB[0][j]][NB[1][j]]<Td)
            Td=T[NB[0][j]][NB[1][j]];
    }
    return(Td);
}
/// Heap_Sort to decrease the time for sorting of the minimum time in narrow-band
void Heap_Sort (struct sv *s)
{
	int i,j,k,dummy,NNP;
	int Ii,Ji,Ik,Jk;

    NNP=*(s->NBP);

	for(j=NNP-2;j>=0;j--)
	{
		k=j;
		while(k<NNP-1)
		{
			i=NNP-1-k;
			i=NNP-1-((i+1)/2-1);
			Ji=*(*(s->NB+1)+i);
			Ii=*(*(s->NB)+i);
			Jk=*(*(s->NB+1)+k);
			Ik=*(*(s->NB)+k);

			if(*(*(s->T+Ji)+Ii)>*(*(s->T+Jk)+Ik))
			{
				dummy=*(*(s->I+Ji)+Ii);
				*(*(s->I+Ji)+Ii) = *(*(s->I+Jk)+Ik);
				*(*(s->I+Jk)+Ik) = dummy;

				dummy=*(*(s->NB)+i);
				*(*(s->NB)+i)=*(*(s->NB)+k);
				*(*(s->NB)+k)=dummy;

				dummy=*(*(s->NB+1)+i);
				*(*(s->NB+1)+i)=*(*(s->NB+1)+k);
				*(*(s->NB+1)+k)=dummy;

				k=i;
			}
			else
				k=NNP-1;
		}
	}

}

int recompute(struct mv *m,struct sv *s,int ip,int jp,int NGx,int NGy,float Dx,float Dy)
{
	int N,i,indc;
	int schemeF(struct mv *m,struct sv *s,int ip,int jp,int xn,int yn,
            float Dx,float Dy);

	indc=schemeF(m,s,ip,jp,NGx,NGy,Dx,Dy);

    if(indc==1)
    {
        if(*(*(s->I+jp)+ip)==Far_pt)
        {
            *(s->NBP)=*(s->NBP)+1;

//Reallocate memory for NB
    		for(i=0;i<2;i++)
        		*(s->NB+i)=(int *) realloc (*(s->NB+i),(*(s->NBP))*2*sizeof(int));

            N = *(s->NBP)-1;
            *(*(s->NB)+N)=ip;
            *(*(s->NB+1)+N)=jp;
            *(*(s->I+jp)+ip)=jp;
            *(*(s->NIt+jp)+ip)=*(*(s->NIt+jp)+ip) + 1;
        }
        else
            *(*(s->NIt+jp)+ip)=*(*(s->NIt+jp)+ip) + 1;
        return(1);
    }
    else
        return(0);
}

// solution solves the quadratic equation for T. If the discriminant becomes negative then it exits
double solution(float F,float vx,float vy, double a, double b, double c,
        double d,int i,int j)
{
	double v[2],FF;
	double Sol[2];
	double dummy1,dummy2,T;
	double coef1,coef2,coef3,discr;

    v[0]=vx;v[1]=vy;
	coef1=0.0;coef2=0.0;coef3=0.0;
    FF=F*F;

	dummy1=FF-v[0]*v[0];
	dummy2=FF-v[1]*v[1];

	coef1=dummy1*a*a-2.0*v[0]*v[1]*a*c+dummy2*c*c;
	coef2=2.0*(dummy1*a*b+v[0]*a-v[0]*v[1]*(a*d+b*c)+v[1]*c+dummy2*c*d);
	coef3=dummy1*b*b+2.0*v[0]*b-2.0*v[0]*v[1]*b*d+2.0*v[1]*d+dummy2*d*d-1.0;

    discr=coef2*coef2-4.0*coef1*coef3;
	if(discr<0)
	{
//		printf("Error: Solution becomes complex at (i,j)= (%d,%d). \n",i,j);
 //       exit(0);
        return(Infinity);
	}

	if(coef1==0)
    {
        printf("coef1 = 0 at (%d, %d) with dummy1 = %lf, dummy2 = %lf  a=%lf c=%lf.\n",i,j,dummy1,dummy2,a,c);
        exit(0);
        return(Infinity);
    }

    discr=sqrt(discr);
    T=2.0*coef1;

    Sol[0]=(-1*coef2+discr)/T;
	Sol[1]=(-1*coef2-discr)/T;
//-------------------------------------------------------------------------------------------
//--------------------------------Choosing the maximum---------------------------------------
//-------------------------------------------------------------------------------------------

	T=((Sol[0]>Sol[1])?Sol[0]:Sol[1]);

	return(T);
} // end of solution

float CharDir(float F, float vx, float vy, double T,double a,double b,double c,double d)
{
	float phi,DTx,DTy;
	float N,D,G;
    double FF;

	DTx=a*T+b;
    DTy=c*T+d;

    FF=F*F;
    G = (1.0-vx*DTx-vy*DTy);

    N = FF*DTy + vy*G;
    D = FF*DTx + vx*G;
    phi = atan2(N,D);

	return(phi);

}// end of function CharDir


// 'Scheme_8_First_Flow'computes the value of T at (i,j) and returns the characteristic direction at (i,j)
int schemeF(struct mv *m,struct sv *s,int i,int j,int NGx,int NGy,float Dx,float Dy)
{

    double Txd,Tyd;

    float F,vx,vy,Norm;
    float e[6][2],Ediff;

	double Tdummy,Tx[4],Ty[4];
    double Ta[6],Tb[6],Tc[6],Td[6];

    double TV1,TV2,Tval[6],TvalCh[6];
    float Phi,Ch1,Ch2;
    double a,b,c,d;
    float M[6][2][2];
	int ip,im,jp,jm,k;
    int Tind;

//The following are TriangleIndicator, TriangleCount,
//TriangleCharacteristicIndicator, TriangleCharacteristicCount.
    int TrI[6],TrC=0,TrChI[6],TrChC=0,TrCComplex=0;

    float CharDir(float F, float vx, float vy, double T,double a,double b,double c,double d);
	double solution(float F,float vx,float vy, double a, double b, double c,
            double d,int i,int j);

	F=*(*(m->F+j)+i);vx=*(*(m->vx+j)+i);vy=*(*(m->vy+j)+i);
    Norm=sqrt(Dx*Dx+Dy*Dy);

    ip=i+1;im=i-1;
    jp=j+1;jm=j-1;
//------------------------------------
// Triangle Coefficients
// Triangle 1:
    if((*(*(s->I+j)+ip)!=Far_pt) && (*(*(s->I+jp)+i)!=Far_pt))
    {
//         M[0][0][0]=-1;M[0][0][1]=0;
//         M[0][1][0]=0;M[0][1][1]=-1;

        TV1=*(*(s->T+j)+ip);
        TV2=*(*(s->T+jp)+i);
        a=-1/Dx;b=TV1/Dx;
        c=-1/Dy;d=TV2/Dy;
  //      printf("T1 (%d %d\n",i,j);
        Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
        if(Tdummy!=Infinity)
        {
            Tval[TrC]=Tdummy;
            Phi=CharDir(F,vx,vy,Tval[TrC],a,b,c,d);
            Ch1=cos(Phi);
            Ch2=sin(Phi);
            if(Ch1<=0 && Ch2<=0)
            {
                TrChI[TrChC]=0;
                TvalCh[TrChC]=Tval[TrC];
                TrChC++;
            }
    //        printf("T1\n");
            TrI[TrC]=0;
            TrC++;
        }
        else
        {
            TrI[TrCComplex]=-1;
            TrCComplex++;
        }

    }
// Triangle 2:

    if((*(*(s->I+jp)+i)!=Far_pt) && (*(*(s->I+jp)+im)!=Far_pt) )
    {
//         M[1][0][0]=0;M[1][0][1]=-1;
//         M[1][1][0]=Dx/Norm;M[1][1][1]=-1*Dy/Norm;

        TV1=*(*(s->T+jp)+i);
        TV2=*(*(s->T+jp)+im);
        a=0;b=(TV1-TV2)/Dx;
        c=-1/Dy;d=TV1/Dy;
 //       printf("T2 (%d %d\n",i,j);
        Tdummy=solution(F,vx,vy,a,b,c,d,i,j);

        if(Tdummy!=Infinity)
        {
            Tval[TrC]=Tdummy;
            Phi=CharDir(F,vx,vy,Tval[TrC],a,b,c,d);
            Ch1=cos(Phi);Ch2=sin(Phi);
            if(Ch1>=0 && Ch2<=(-1*Dy*Ch1/Dx))
            {
                TrChI[TrChC]=1;
                TvalCh[TrChC]=Tval[TrC];
                TrChC++;
            }
            TrI[TrC]=1;TrC++;

   //         printf("T2\n");
        }
        else
        {
            TrI[TrCComplex]=-2;
            TrCComplex++;
        }
                           // printf("Hello %d %d  2\n",i,j);
    }
// Triangle 3:

    if((*(*(s->I+jp)+im)!=Far_pt) && (*(*(s->I+j)+im)!=Far_pt))
    {
//         M[2][0][0]=Dx/Norm;M[2][0][1]=-1*Dy/Norm;
//         M[2][1][0]=1;M[2][1][1]=0;

        TV1=*(*(s->T+jp)+im);
        TV2=*(*(s->T+j)+im);
        a=1/Dx;b=-TV2/Dx;
        c=0;d=(TV1-TV2)/Dy;
   //     printf("T3 (%d %d\n",i,j);
        Tdummy=solution(F,vx,vy,a,b,c,d,i,j);

        if(Tdummy!=Infinity)
        {
            Tval[TrC]=Tdummy;
            Phi=CharDir(F,vx,vy,Tval[TrC],a,b,c,d);
            Ch1=cos(Phi);Ch2=sin(Phi);
            if(Ch1>=(-1*Dx*Ch2/Dy) && Ch2<=0)
            {
                TrChI[TrChC]=2;
                TvalCh[TrChC]=Tval[TrC];
                TrChC++;
            }

            TrI[TrC]=2;TrC++;
   //         printf("T3\n");
        }
        else
        {
            TrI[TrCComplex]=-3;
            TrCComplex++;
        }
    }
// Triangle 4:

    if((*(*(s->I+j)+im)!=Far_pt) && (*(*(s->I+jm)+i)!=Far_pt))
    {
//         M[3][0][0]=1;M[3][0][1]=0;
//         M[3][1][0]=0;M[3][1][1]=1;

        TV1=*(*(s->T+j)+im);
        TV2=*(*(s->T+jm)+i);
        a=1/Dx;b=-TV1/Dx;
        c=1/Dy;d=-TV2/Dy;
  //      printf("T4 (%d %d\n",i,j);
        Tdummy=solution(F,vx,vy,a,b,c,d,i,j);

        if(Tdummy!=Infinity)
        {
            Tval[TrC]=Tdummy;
            Phi=CharDir(F,vx,vy,Tval[TrC],a,b,c,d);
            Ch1=cos(Phi);Ch2=sin(Phi);
            if(Ch1>=0 && Ch2>=0)
            {
                TrChI[TrChC]=3;
                TvalCh[TrChC]=Tval[TrC];
                TrChC++;
            }

            TrI[TrC]=3;TrC++;
   //         printf("T4\n");
        }
        else
        {
            TrI[TrCComplex]=-4;
            TrCComplex++;
        }
    }
// Triangle 5:

    if((*(*(s->I+jm)+i)!=Far_pt) && (*(*(s->I+jm)+ip)!=Far_pt))
    {
//         M[4][0][0]=0;M[4][0][1]=1;
//         M[4][1][0]=-1*Dx/Norm;M[4][1][1]=Dy/Norm;

        TV1=*(*(s->T+jm)+i);
        TV2=*(*(s->T+jm)+ip);
        a=0;b=(TV2-TV1)/Dx;
        c=1/Dy;d=-TV1/Dy;
     //   printf("T5 (%d %d\n",i,j);
        Tdummy=solution(F,vx,vy,a,b,c,d,i,j);

        if(Tdummy!=Infinity)
        {
            Tval[TrC]=Tdummy;
            Phi=CharDir(F,vx,vy,Tval[TrC],a,b,c,d);
            Ch1=cos(Phi);Ch2=sin(Phi);
            if(Ch1<=0 && Ch2>=(-1*Dy*Ch1/Dx))
            {
                TrChI[TrChC]=4;
                TvalCh[TrChC]=Tval[TrC];
                TrChC++;
            }

            TrI[TrC]=4;TrC++;
  //          printf("T5\n");
        }
        else
        {
            TrI[TrCComplex]=-5;
            TrCComplex++;
        }
    }
// Triangle 6:

    if((*(*(s->I+jm)+ip)!=Far_pt) && (*(*(s->I+j)+ip)!=Far_pt))
    {
//         M[5][0][0]=-1*Dx/Norm;M[5][0][1]=Dy/Norm;
//         M[5][1][0]=-1;M[5][1][1]=0;

        TV1=*(*(s->T+jm)+ip);
        TV2=*(*(s->T+j)+ip);
        a=-1/Dx;b=TV2/Dx;
        c=0;d=(TV2-TV1)/Dy;
  //      printf("T6 (%d %d\n",i,j);
        Tdummy=solution(F,vx,vy,a,b,c,d,i,j);

        if(Tdummy!=Infinity)
        {
            Tval[TrC]=Tdummy;
            Phi=CharDir(F,vx,vy,Tval[TrC],a,b,c,d);
            Ch1=cos(Phi);Ch2=sin(Phi);
            if(Ch1<=(-1*Dx*Ch2/Dy) && Ch2>=0)
            {
                TrChI[TrChC]=5;
                TvalCh[TrChC]=Tval[TrC];
                TrChC++;
            }

            TrI[TrC]=5;
            TrC++;
    //        printf("T6\n");
        }
        else
        {
            TrI[TrCComplex]=-6;
            TrCComplex++;
        }
    }
//------------------------------------

/* At this stage the available triangles are Computed.
 * Now we have to look for the Difficulties.
 */

//Dificulty 1:

    if(TrC==0 && TrCComplex==0)
    {
//Triangle 1: Case 1
        if(*(*(s->I+j)+ip)!=Far_pt && (*(*(s->I+jp)+i)==Far_pt))
        {
            TV1=*(*(s->T+j)+ip);

            a=-1/Dx;b=TV1/Dx;
            c=0;d=0;

            Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
            if(Tdummy!=Infinity)
            {
                Tval[TrC]=Tdummy;
                TrI[TrC]=0;TrC++;
                printf("T1 E1\n");
            }
            else
            {
                TrI[TrCComplex]=-1;
                TrCComplex++;
            }
        }
//Triangle 1: Case 2
        else if (*(*(s->I+j)+ip)==Far_pt && (*(*(s->I+jp)+i)!=Far_pt))
        {
            TV2=*(*(s->T+jp)+i);
            a=0;b=0;
            c=-1/Dy;d=TV2/Dy;

            Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
            if(Tdummy!=Infinity)
            {
                Tval[TrC]=Tdummy;
                TrI[TrC]=0;TrC++;
                printf("T1 E2\n");
            }
            else
            {
                TrI[TrCComplex]=-1;
                TrCComplex++;
            }
        }
//Triangle 2: Case 1
        if ((*(*(s->I+jp)+i)!=Far_pt) && *(*(s->I+jp)+im)==Far_pt)
        {
            TV1=*(*(s->T+j)+ip);
            a=-1/Dx;b=TV1/Dx;
            c=-1/Dy;d=TV1/Dy;
            Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
            if(Tdummy!=Infinity)
            {
                Tval[TrC]=Tdummy;
                TrI[TrC]=1;TrC++;
                printf("T2 E2\n");
            }
            else
            {
                TrI[TrCComplex]=-2;
                TrCComplex++;
            }
        }
//Triangle 2: Case 2
        else if ((*(*(s->I+jp)+i)==Far_pt) && *(*(s->I+jp)+im)!=Far_pt)
        {
            TV2=*(*(s->T+jp)+im);
            a=1/Dx;b=-TV2/Dx;
            c=0;d=0;
            Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
            if(Tdummy!=Infinity)
            {
                Tval[TrC]=Tdummy;
                TrI[TrC]=1;TrC++;
                printf("T2 E3\n");
            }
            else
            {
                TrI[TrCComplex]=-2;
                TrCComplex++;
            }
        }
//Triangle 3: Case 1
        if (*(*(s->I+jp)+im)!=Far_pt && *(*(s->I+j)+im)==Far_pt)
        {
            TV1=*(*(s->T+jp)+im);
            a=0;b=0;
            c=-1/Dy;d=TV1/Dy;
            Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
            if(Tdummy!=Infinity)
            {
                Tval[TrC]=Tdummy;
                TrI[TrC]=3;TrC++;
                printf("T3 E3\n");
            }
            else
            {
                TrI[TrCComplex]=-3;
                TrCComplex++;
            }
        }
//Triangle 3: Case 2
        else if (*(*(s->I+jp)+im)==Far_pt && *(*(s->I+j)+im)!=Far_pt)
        {
             TV2=*(*(s->T+j)+im);
             a=1/Dx;b=-TV2/Dx;
             c=1/Dy;d=-TV2/Dy;

             Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
             if(Tdummy!=Infinity)
             {
                Tval[TrC]=Tdummy;
                TrI[TrC]=2;TrC++;
                printf("T3 E4\n");
             }
             else
             {
                 TrI[TrCComplex]=-3;
                 TrCComplex++;
             }
        }
//Triangle 4: Case 1
        if (*(*(s->I+j)+im)!=Far_pt && *(*(s->I+jm)+i)==Far_pt)
        {
            TV1=*(*(s->T+j)+im);

            a=1/Dx;b=-TV1/Dx;
            c=0;d=0;

            Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
            if(Tdummy!=Infinity)
            {
                Tval[TrC]=Tdummy;

                TrI[TrC]=4;TrC++;
                printf("T4 E4\n");
            }
            else
            {
                TrI[TrCComplex]=-4;
                TrCComplex++;
            }
        }
//Triangle 4: Case 2
        else if (*(*(s->I+j)+im)==Far_pt && *(*(s->I+jm)+i)!=Far_pt)
        {
            TV2=*(*(s->T+jm)+i);

            a=0;b=0;
            c=1/Dy;d=-TV2/Dy;

            Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
            if(Tdummy!=Infinity)
            {
                Tval[TrC]=Tdummy;

                TrI[TrC]=4;TrC++;

            }
            else
            {
                TrI[TrCComplex]=-4;
                TrCComplex++;
            }
        }
//Triangle 5: Case 1
        if (*(*(s->I+jm)+i)!=Far_pt && *(*(s->I+jm)+ip)==Far_pt)
        {
            TV1=*(*(s->T+jm)+i);

            a=1/Dx;b=-TV1/Dx;
            c=1/Dy;d=-TV1/Dy;

            Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
            if(Tdummy!=Infinity)
            {
                Tval[TrC]=Tdummy;

                TrI[TrC]=4;TrC++;
                printf("T5 E5\n");
            }
            else
            {
                TrI[TrCComplex]=-5;
                TrCComplex++;
            }
        }
//Triangle 5: Case 2
        else if (*(*(s->I+jm)+i)==Far_pt && *(*(s->I+jm)+ip)!=Far_pt)
        {
            TV2=*(*(s->T+jm)+ip);

            a=-1/Dx;b=TV2/Dx;
            c=0;d=0;

            Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
            if(Tdummy!=Infinity)
            {
                Tval[TrC]=Tdummy;

                TrI[TrC]=4;TrC++;
                printf("T5 E6\n");
            }
            else
            {
                TrI[TrCComplex]=-5;
                TrCComplex++;
            }
        }
//Triangle 6: Case 1
        if (*(*(s->I+jm)+ip)!=Far_pt && *(*(s->I+j)+ip)==Far_pt)
        {
            TV1=*(*(s->T+jm)+ip);

            a=0;b=0;
            c=1/Dy;d=-TV1/Dy;

            Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
            if(Tdummy!=Infinity)
            {
                Tval[TrC]=Tdummy;

                TrI[TrC]=4;TrC++;
                printf("T6 E6\n");
            }
            else
            {
                TrI[TrCComplex]=-6;
                TrCComplex++;
            }
        }
//Triangle 6: Case 2
        else if (*(*(s->I+jm)+ip)==Far_pt && *(*(s->I+j)+ip)!=Far_pt)
        {
            TV2=*(*(s->T+j)+ip);

            a=-1/Dx;b=TV2/Dx;
            c=-1/Dy;d=TV2/Dy;

            Tdummy=solution(F,vx,vy,a,b,c,d,i,j);
            if(Tdummy!=Infinity)
            {
                Tval[TrC]=Tdummy;

                TrI[TrC]=4;TrC++;
                printf("T6 E1\n");
            }
            else
            {
                TrI[TrCComplex]=-6;
                TrCComplex++;
            }
        }
        if(TrC!=0)
        {
            Tdummy=Tval[0];
            for(k=1;k<TrC;k++)
                Tdummy=((Tdummy>Tval[k])?Tdummy:Tval[k]);
            *(*(s->T+j)+i)=Tdummy;
            *(*(s->Q+j)+i)=-1;
            return(1);
        }
        else
        {
            printf("Difficulty 1: Status retained at (%d %d)\n",i,j);
            *(*(s->Q+j)+i)=-1;
            return(0);
        }
    }
    else if (TrC!=0)// && TrCComplex==0)
    {
//Difficulty 2:
        if(TrChC==0)
        {
            Tdummy=Tval[0];
    //        printf("D2 %d %d\n",i,j);
            for(k=1;k<TrC;k++)
                Tdummy=((Tdummy<Tval[k])?Tdummy:Tval[k]);
            *(*(s->T+j)+i)=Tdummy;
            *(*(s->Q+j)+i)=1;
            return(1);
        }
        else if (TrChC>1)
        {
//Difficulty 3:
    //        printf("D3 %d %d\n",i,j);
            Tdummy=TvalCh[0];
            for(k=1;k<TrChC;k++)
                Tdummy=((Tdummy<TvalCh[k])?Tdummy:TvalCh[k]);
            *(*(s->T+j)+i)=Tdummy;
            *(*(s->Q+j)+i)=1;
            return(1);
        }
//No difficulty, everything is fine.
        else
        {
     //       printf("All right %d %d\n",i,j);
            *(*(s->T+j)+i)=TvalCh[0];//This case is when TrChC=1, therefore 1 less is the index for TvalCh
            *(*(s->Q+j)+i)=1;
            return(1);
        }
    }
    else
    {
        printf("Problem at (%d, %d).  All triangles returned complex value\n",i,j);
        printf("TrC=%d\n",TrC);
        for(k=1;k<6;k++)
            printf("%lf\n",Tval[k]);
        exit(0);
        return(0);
    }

}
