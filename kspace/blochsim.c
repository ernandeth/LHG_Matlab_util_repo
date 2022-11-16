#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	MI	prhs[0]
#define	BEFF	prhs[1]
#define TONE    prhs[2]
#define TTWO    prhs[3]
#define DT      prhs[4]
#define NSTEP   prhs[5]

/* Output Arguments */

#define	MOUT	plhs[0]

#define PI 3.14159265

static void blochsim(double *Mi,
		     double *beff,
		     double T1,
		     double T2,
		     double dt,
		     int nstep,
		     double *M)
{
    
  /* dt in ms */
  double gambar = 42.57e3;
  double gam = gambar*2*PI;
  T1 = dt / T1;
  T2 = (1. - dt / T2);

  /* Put Beff into appropriate units */
  int i = 0;
  for(i=0;i<nstep;i++)
    {
      
      beff[0*nstep+i] = beff[0*nstep+i]*dt*gam;
      beff[1*nstep+i] = beff[1*nstep+i]*dt*gam;
      beff[2*nstep+i] = beff[2*nstep+i]*dt*gam;
      
    }
  
  M[0*(nstep+1)]=Mi[0];
  M[1*(nstep+1)]=Mi[1];
  M[2*(nstep+1)]=Mi[2];
  
  int lp = 0;
  for(lp=1;lp<=nstep;lp++)
    {
      
      double Bx = beff[0*nstep+lp-1]; 
      double By = beff[1*nstep+lp-1]; 
      double Bz = beff[2*nstep+lp-1];

      /*
	Compute sines and cosines of field angles:
	Theta = angle w.r.t. positive z axis
	Phi = angle w.r.t. positive x axis
	Psi = angle w.r.t. transformed positive x axis
      */

      double Bmag = sqrt(Bx*Bx+By*By+Bz*Bz);  /* Magnitude of applied field */
      double Btrans = sqrt(Bx*Bx+By*By);      /* Magnitude of transverse applied field */

      double ct = 1;
      if(Bmag > 0)
	ct = Bz/Bmag;  /* cos(theta) */

      double st = sqrt(1 - ct*ct);  /* sin(theta) > 0 */

      double cphi = 1;
      if(Btrans > 0)
	cphi = Bx/Btrans;  /* cos(phi) */
      
      double sphi; 
      if(By < 0.)
	sphi = sqrt(1 - cphi*cphi)*(-1.);
      else
	sphi = sqrt(1 - cphi*cphi);

      double cpsi = cos(Bmag);  /* cos(psi) */
      double spsi = sin(Bmag);  /* sin(psi) */
      
      double Mx1,My1,Mz1;
      
      if(Bmag > 0)
	{
	  double Mx0 = M[0*(nstep+1)+lp-1];
	  double My0 = M[1*(nstep+1)+lp-1];
	  double Mz0 = M[2*(nstep+1)+lp-1];
	  
	  Mx1 = cphi*(ct*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0))+st*(ct*Mz0+st*(sphi*My0+cphi*Mx0))) - sphi*(-1.*spsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + cpsi*(cphi*My0-sphi*Mx0));
	  My1 = sphi*(ct*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0))+st*(ct*Mz0+st*(sphi*My0+cphi*Mx0))) + cphi*(-1.*spsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + cpsi*(cphi*My0-sphi*Mx0));
	  Mz1 = ct*(ct*Mz0+st*(sphi*My0+cphi*Mx0)) - st*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0));
	}
      else 
	{
	  
	  Mx1 = M[0*(nstep+1)+lp-1];
	  My1 = M[1*(nstep+1)+lp-1];
	  Mz1 = M[2*(nstep+1)+lp-1];
	  
	}
      
      /* relaxation effects: "1" in Mz since Mo=1 by assumption */
      
      M[0*(nstep+1)+lp] = Mx1*T2;
      M[1*(nstep+1)+lp] = My1*T2;
      M[2*(nstep+1)+lp] = Mz1 + (1 - Mz1)*T1;
      
    }
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *Mi; 
    double *beff;
    double *M;

    /* Check for proper number of arguments */
    
    if (nrhs != 6) { 
	mexErrMsgTxt("Six input arguments required."); 
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 

    /* Assign pointers to the various parameters */ 
    Mi = mxGetPr(MI);
    beff = mxGetPr(BEFF); 

    double T1 = mxGetScalar(TONE);
    double T2 = mxGetScalar(TTWO);
    double dt = mxGetScalar(DT);
    int nstep = (int)mxGetScalar(NSTEP);
    
    /* Create a matrix for the return argument */ 
    MOUT = mxCreateDoubleMatrix(nstep+1, 3, mxREAL); 
    M = mxGetPr(MOUT);
    
    /* Do the actual computations in a subroutine */
    blochsim(Mi,beff,T1,T2,dt,nstep,M);
    return;
    
}
