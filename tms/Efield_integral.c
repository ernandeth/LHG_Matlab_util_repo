/*
 *  Efield_integral.c
 *  
 *
 *  Created by luis hernandez-garcia on Tue May 31 2005.
 *  Copyright (c) 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "mex.h"
    /*Permeability H/m*/


void mexFunction(int nlhs,
                mxArray *plhs[],
                int nrhs,
                const mxArray *prhs[])
{


	/* Check for proper number of arguments */
    
	if (nrhs != 9) 
	{ 
		mexErrMsgTxt("\nUSAGE: Efield_integral(Ex, Ey,Ez,wire,dwire,dI_dt,x,y,z)"); 
	}
	else if (nlhs > 0) 
	{
		mexErrMsgTxt("Too many output arguments."); 
	} 
	
	#define	MU  1.26e-6
	#define	PI 3.14159

    /* inputs from matlab*/
    double *pdEx = mxGetPr(prhs[0]);
    double *pdEy = mxGetPr(prhs[1]);
    double *pdEz = mxGetPr(prhs[2]);
    double *pdWire = mxGetPr(prhs[3]);
    double *pdDwire = mxGetPr(prhs[4]);
    double dDI_dt = mxGetScalar(prhs[5]);
    double *pdX = mxGetPr(prhs[6]);
    double *pdY = mxGetPr(prhs[7]);
    double *pdZ = mxGetPr(prhs[8]);
    
    /* No outputs to matlab - we'll work on the pointers
    double *result;
    plhs[0] = mxCreateDoubleMatrix(tLen, 1, mxREAL);
    result=mxGetPr(plhs[0]);
    */

    
    int iVox;
    int	iVOXELS = mxGetM(prhs[0]);
    double dR;
    int	iNsegs = mxGetM(prhs[3]);
    int iSeg;
    double dK ;

	dK =  - dDI_dt * MU / (4 * PI);

    
    for (iVox=0; iVox<=iVOXELS; iVox++)
    {
        pdEx[iVox] = 0.0;
        pdEy[iVox] = 0.0;
        pdEz[iVox] = 0.0;

        // intetgrate over the length of the wire
        for (iSeg=0; iSeg<=iNsegs; iSeg++)
        {
            // calculate the distance from voxel to each wire segment
            dR = sqrt( 
               (pdX[iVox] - pdWire[0*iNsegs + iSeg]) * (pdX[iVox] - pdWire[0*iNsegs + iSeg]) 
             + (pdY[iVox] - pdWire[1*iNsegs + iSeg]) * (pdY[iVox] - pdWire[1*iNsegs + iSeg]) 
             + (pdZ[iVox] - pdWire[2*iNsegs + iSeg]) * (pdZ[iVox] - pdWire[2*iNsegs + iSeg])); 
             
            // E field  = integral{ dI_dt(l) / distance(l)  }  
            pdEx[iVox] += dK * pdDwire[0*iNsegs + iSeg] / dR;
            pdEy[iVox] += dK * pdDwire[1*iNsegs + iSeg] / dR;
            pdEz[iVox] += dK * pdDwire[2*iNsegs + iSeg] / dR;
        }
        
    
    }

}
