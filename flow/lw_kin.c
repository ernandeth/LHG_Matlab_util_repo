#include <math.h>
#include "mex.h"

void mexFunction( 	int nlhs, 
			mxArray *plhs[], 
		  	int nrhs, 
			const mxArray *prhs[] )
{


/* Check for proper number of arguments */
    
if (nrhs != 12) 
{ 
	mexErrMsgTxt("\nUSAGE: lw_kin(tvec, f,V,A,tis,art,lambda,R1a,R1t,dist)"); 
} 
else if (nlhs > 1) 
{
	mexErrMsgTxt("Too many output arguments."); 
} 


int t, x,i;
int tLen, xLen, nrows, ncols;
int xBin;
double dx, dt;
double Vt, ft, Ax, At, Axx, Att, artt,  Tt, Ttt;
double *tmp;

/* inputs from matlab*/
double *tvec = mxGetPr(prhs[0]);
double *f = mxGetPr(prhs[1]);
double *V = mxGetPr(prhs[2]);
double *A = mxGetPr(prhs[3]);
double *tis = mxGetPr(prhs[4]); 
double *art = mxGetPr(prhs[5]);
double lambda = mxGetScalar(prhs[6]);
double R1a = mxGetScalar(prhs[7]);
double R1t = mxGetScalar(prhs[8]);
double dist = (double)mxGetScalar(prhs[9]);
int    CM = (int)mxGetScalar(prhs[10]);
int    SECONDS = (int)mxGetScalar(prhs[11]);
/*
double *dummy = mxGetPr(prhs[10]);
*/

tLen = mxGetN(prhs[0]);
xLen = mxGetM(prhs[3]);
/*
printf("%\n tLen ...%d time points..\n",tLen);
printf("%\n xLen ...%d space points..\n",xLen);
*/
/* outputs to matlab :
double *result;
plhs[0] = mxCreateDoubleMatrix(tLen, 1, mxREAL);
result=mxGetPr(plhs[0]);
*/
nrows=xLen;
ncols=tLen;

dx=1.0/CM;
dt=1.0/SECONDS;

xBin = (int)(round(CM*dist));
/*
mxArray *rhs[1];
rhs[0] = mxCreateDoubleMatrix(xLen,tLen,mxREAL);
tmp=mxGetPr(rhs[0]);

for (i=0;i<tLen;i++) tmp[i]=tis[i];
mexCallMATLAB(0,NULL,1,rhs,"plot");
*/

for (t=1; t < tLen-1; t++)
{
    /* first time derivative of V and f:*/
    Vt = ( V[t+1] - V[t-1] ) / (2*dt);
    ft = ( f[t+1] - f[t-1] ) / (2*dt);


    for (x=1; x < xLen-1; x++)
    { 
        /*first derivatives in space and time*/
        Ax = ( A[(x+1) + nrows*t ] - A[(x-1)+ nrows*t ] ) / (2*dx);
        At = -V[t]*Ax - R1a*A[x + t*nrows ];
        /*second derivative in space and time*/
        Axx = (A[(x-1) + nrows*t ] - 2*A[x + nrows*t] + A[(x+1) + nrows*t]) / (dx*dx);
        Att = -V[t]*(-V[t]*Axx - R1a*Ax) - Vt*Ax - R1a*At;

        /*compute the next step (Taylor):*/
        A[x + nrows*(t+1)] = A[x + nrows*t] + dt*At + (dt*dt)*Att/2 ;

        if (isnan(A[x + nrows*t]))
        {
             printf("\nA=nan  at t=%d	x=%d", t,x);
             printf("\nt= %d  art= %f  f=%f  V=%f   Tt= %f ft= %f  tis= %f", 
                   t, art[t], f[t], V[t], Tt, ft,tis[t]);
             return;
         }
  
    }
/*
    for (i=0;i<xLen*tLen;i++) tmp[i]=A[i];
    mexCallMATLAB(0,NULL,1,rhs,"imagesc");
    mexCallMATLAB(0,NULL,0,NULL,"drawnow");
*/
    /*the tissue compartment is only this one element at "dist" cm.*/
    art[t] =  A[xBin + nrows*t]  - f[t] * A[xBin + nrows*t ]  ;
    A[xBin + nrows*t] = A[xBin + nrows*t]  - f[t] * A[xBin + nrows*t ]  ;
    
    
    /* the tissue compartment using same method....*/
    artt = (art[t]-art[t-1])/dt; 
    Tt = f[t]*art[t] - R1t*tis[t];
    Ttt = ft*art[t] + f[t]*artt - R1t*Tt;
    
    tis[t+1] = tis[t] + dt*Tt +(dt*dt)*Ttt/2;

    if (isnan(tis[t]))
    {
         printf("\ntis=nan at t=%d	x=%d", t,x);
         printf("\nt= %d  art= %f  f=%f Tt= %f ft= %f  tis= %f", 
	        t, art[t], f[t], Tt, ft,tis[t]);
    }
  
}

/*
for(i=0;i<100;i++) dummy[i]=3*dummy[i];
*/

/*return;*/
}
