/* charfmmnv.c -> CHARacteristic Fast Marching Method New Version
 *
 *What is new in this?
 * 
 * 1. The elaborate coding of the previous version has been greatly reduced.
 * 2. Upgraded to Structures.
 * 3. Input taken from an external input file named 'input.dat'.
 * 4. All the data files (both input and output) are kept in a seperate
 *    directory called 'Data'.
 */ 
//---------------Standard C librabries are included here-------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "mex.h"
//-------------------------------------------------------------------------------------
//-------This is a user defined header file available in the same directory------------

#include "Tcharfmmnv.h"
//-------------------------------------------------------------------------------------
//void mexFunction( int nlhs, mxArray *plhs[], 
//		  int nrhs, const mxArray *prhs[] ) 
//{
int main () {
//---------------------Variables Definition--------------------------------------------
	void charfmm_without_ref(int NGx, int NGy, int NAP, float Dx, float Dy,
        struct mv *m,struct sv *s);
    
	double xmin, xmax, ymin, ymax, npoints,vy, vx, inds,indw;
	int npts,ITs,ITw;
    
    int i,j,k,check;
    int NGx,NGy,NAP,NNP;
    double Tv;
    float Dx,Dy;
    
    struct sv *s;
	struct mv *m;
    
    FILE *fpi;
//-------------------------------------------------------------------------------------
    fpi=fopen("Data/inputN.dat","r");
    fscanf(fpi,"%d%d%d%d",&NGx,&NGy,&NAP,&NNP);
    fscanf(fpi,"%f%f",&Dx,&Dy);
    fclose(fpi);
//----------------------------Memory Allocation----------------------------------------    
    m = (struct mv *) malloc (sizeof(struct mv));
	s = (struct sv *) malloc (sizeof(struct sv));
	check=mmalloc(m,NGx,NGy);
	if(check==0)
	{
		printf("Insufficent memory. Program terminated\n");
        exit(0);
	}
	check=smalloc(s,NGx,NGy,NNP);
	if(check==0)
	{
		printf("Insufficent memory. Program terminated\n");
        exit(0);
	}
    *(s->NBP)=NNP;
//-------------------------------------------------------------------------------------
//------------------------------- Assigning F and v -----------------------------------
    fpi=fopen("Data/inputD.dat","r");
    for(j=1;j<=NGy;j++)
	{
		for(i=1;i<=NGx;i++)
			fscanf(fpi,"%f%f%f",*(m->F+j)+i,*(m->vx+j)+i,*(m->vy+j)+i);
	}
//-------------------------------------------------------------------------------------
//-------------------Initialization of Alive Points and Narrow Band Points-------------
//--------------------------------[ALIVE POINTS]---------------------------------------
	for(k=0;k<NAP;k++)
	{
		fscanf(fpi,"%d%d%lf",&i,&j,&Tv);
		*(*(s->T+j)+i)=Tv;
		*(*(s->I+j)+i)=Alive_pt;
	}
//-------------------------------------------------------------------------------------
//---------------------------------[NARROW BAND POINTS]--------------------------------
    for(k=0;k<NNP;k++)
	{
		fscanf(fpi,"%d%d%lf",&i,&j,&Tv);
		*(*(s->T+j)+i)=Tv;
		*(*(s->I+j)+i)=j;
		*(*(s->NB)+k)=i;        //Index of the NB array:  0 --> x position.
		*(*(s->NB+1)+k)=j;      //Index of the NB array:  1 --> y position.
	}
//-------------------------------------------------------------------------------------
    fclose(fpi);
    
//-------------------------------------------------------------------------------------
    charfmm_without_ref(NGx, NGy, NAP, Dx, Dy, m, s);
	return(0);
}

//-------------------------Main program starts here------------------------------------
void charfmm_without_ref(int NGx, int NGy, int NAP, float Dx, float Dy,
        struct mv *m,struct sv *s)
{
//---------------------------Variables are declared here-------------------------------
	int i,j,I,J,NNP,indc;
    FILE *fp1, *fp2, *fp3, *fp;


//-------------------------------------------------------------------------------------
//|------------------------------MARCHING FORWARD-------------------------------------|
//=====================================================================================
	while(NAP<=NGx*NGy-1) 
	{
        NNP=*(s->NBP);
//  printf("NAP = %d\n",NAP);   
        if(NNP==0)     //If everything is fine, this will never happen.
        {
            printf("WARNING: Unexpected Termination of the computational loop.\n");
            break;              
        }
//-------------------------------------------------------------------------------------
//STEP: 1  : Use Heap Sort to obtain the minimum value among Narrow Band Points
        Heap_Sort(s);
		J=*(*(s->NB+1)+(NNP-1));
		I=*(*(s->NB)+(NNP-1));
     
//STEP: 2  : Marking Alive Points from the Narrow Band set
 //To get unstable solution comment the dashed portion
//------------------------------------------------------------------------------------ 
       
        if(*(*(s->Q+J)+I)<0)
        {
            
            indc=schemeF(m,s,I,J,NGx,NGy,Dx,Dy);  
            if(indc==1)
            {
                NAP++;
                *(*(s->I+J)+I)=Alive_pt;
                *(s->NBP)=*(s->NBP)-1;
            }
        }
 
//------------------------------------------------------------------------------------
        else
        {
            NAP++;
            *(*(s->I+J)+I)=Alive_pt;
            *(s->NBP)=*(s->NBP)-1;
        }
               
//STEP: 3  :
//-------------------------------------------------------------------------------------
//------------------------Calculating the solution at Neighbor Points------------------
//-------------------------------------------------------------------------------------

        if(I<NGx) 
        {   
            i=I+1;j=J;  
            if (*(*(s->I+j)+i)!=Alive_pt)
                recompute(m,s,i,j,NGx,NGy,Dx,Dy);
        }
        if(J<NGy) 
        {   
            i=I;j=J+1;  
            if (*(*(s->I+j)+i)!=Alive_pt)
                recompute(m,s,i,j,NGx,NGy,Dx,Dy);
        }
        
         if(I>1 && J<NGy) 
        {   
            i=I-1;j=J+1;  
            if (*(*(s->I+j)+i)!=Alive_pt)
                recompute(m,s,i,j,NGx,NGy,Dx,Dy);
        }
        if(I>1) 
        {   
            i=I-1;j=J;  
            if (*(*(s->I+j)+i)!=Alive_pt)
                recompute(m,s,i,j,NGx,NGy,Dx,Dy);
        }
        if(J>1) 
        {   
            i=I;j=J-1;  
            if (*(*(s->I+j)+i)!=Alive_pt)
                recompute(m,s,i,j,NGx,NGy,Dx,Dy);
        }
        
        if(I<NGx && J>1) 
        {   
            i=I+1;j=J-1; 
            if (*(*(s->I+j)+i)!=Alive_pt)
                recompute(m,s,i,j,NGx,NGy,Dx,Dy);
        }
		 
		 
//-------------------------------------------------------------------------------------		
	}
 
    for(i=1;i<NGx+1;i++)
    {
        for(j=1;j<NGy+1;j++)
        {
            if(*(*(s->I+j)+i)==Far_pt)
            {
                printf("This point is still Far Away: %d %d %lf\n",i,j,*(*(s->T+j)+i));
            }
        }
    }
//-----------------------------The Calculation is over---------------------------------------
    printf("The Calculation is finished successfully with NNP = %d\n",*(s->NBP));
//	system("date");
//--------------------------Storing the T values in the file T.dat---------------------------
	fp1=fopen("Data/T_1st.dat","w");
	for(j=1;j<=NGy;j++)
	{
		for(i=1;i<=NGx;i++)
		 {	
			fprintf(fp1,"%lf ",*(*(s->T+j)+i));
		 }
		fprintf(fp1,"\n");
	}
	fclose(fp1);

//--------------------------Storing the Q1 values in the file C.dat---------------------------

	fp2=fopen("Data/Q_1st.dat","w");
	for(j=1;j<=NGy;j++)
	{
		for(i=0;i<NGx;i++)
		 {	
			fprintf(fp2,"%d ",*(*(s->Q+j)+i));
		 }
		fprintf(fp2,"\n");
	}
	fclose(fp2);  
//--------------------------Storing the NIt values in the file C.dat---------------------------

	fp2=fopen("Data/NIt_1st.dat","w");
	for(j=1;j<=NGy;j++)
	{
		for(i=0;i<NGx;i++)
		 {	
			fprintf(fp2,"%d ",*(*(s->NIt+j)+i));
		 }
		fprintf(fp2,"\n");
	}
	fclose(fp2);  


//---------------------------------Program finished------------------------------------
}