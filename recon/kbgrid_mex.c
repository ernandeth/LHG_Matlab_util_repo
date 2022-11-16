/*
 * kbgrid_mex.c
 * Copyright 2005, Douglas C. Noll, University of Michigan
 */
#include "mex.h"
#include "math.h"
#include "stdio.h"
#define mexFail(str) {mexEvalString("sprintf('fail: " #str "')"); return;}
#define mexFail0(str) {mexEvalString("sprintf('fail: " #str "')"); return 0;}

static void kbgrid(
		const double *p_kx,
		const double *p_ky,
		const int ndat,
		const double *r_dat,
		const double *i_dat,
		const int nim,
		const double W,
		const double *p_kernel,
		const int nk,
		double *r_fm,
		double *i_fm)
{
	const double dwin = W/2;
	const double nwby2 = (nk-1)/2;
	int i,lx,ly;
	int kptr;
	double dkx, dky, wx, w;

	for (i = 0; i < ndat; i++ ) {
		dkx = *(p_kx+i) +nim/2;
		dky = *(p_ky+i) +nim/2;
		/* printf("%f %f %d\n",dkx,dky,nim);  */
		for (lx = ceil(dkx-dwin); lx<=floor(dkx+dwin); lx++) {
			if ((lx<0) || (lx>=nim)) continue;
			kptr = (int) rint((((dkx-lx)/dwin)*nwby2 )) + nwby2;
			wx = *(p_kernel + kptr);
			/* printf("%f %d %d %d\n",wx,i,kptr,lx);  */
			for (ly = ceil(dky-dwin); ly<=floor(dky+dwin); ly++) {
				if ((ly<0) || (ly>=nim)) continue;
				kptr = (int) rint((((dky-ly)/dwin)*nwby2 )) + nwby2;
				w = wx * *(p_kernel + kptr);
				/* printf("%f %f %d %d %d\n",w,wx,i,kptr,lx*nim+ly);
				   printf("%f %f\n",*(r_dat + i),*(i_dat + i));
				   printf("%f %f\n",*(r_fm + lx*nim + ly),*(i_fm + lx*nim + ly)); */
				*(r_fm + lx*nim + ly) += *(r_dat + i) * w;
				*(i_fm + lx*nim + ly) += *(i_dat + i) * w;
			}
		}
	}
}


int kbgrid_mex(
		mxArray *plhs[],
		const mxArray *mx_kx,
		const mxArray *mx_ky,
		const mxArray *mx_datr,
		const mxArray *mx_dati,
		const mxArray *mx_nim,
		const mxArray *mx_W,
		const mxArray *mx_kernel)
{
	int i,nn,W,ndat;
	const int ndat1 = mxGetN(mx_kx);
	const int ndat2 = mxGetM(mx_kx);
	const int ndim = mxGetNumberOfDimensions(mx_kx);
	const int N = (ndim > 2) ?  (mxGetDimensions(mx_kx))[2] : 1;
	const int nk = mxGetM(mx_kernel);	/* # of kernel samples */

	const double *nim = mxGetPr(mx_nim);
	const double *www = mxGetPr(mx_W);

	const double *p_kx = mxGetPr(mx_kx);
	const double *p_ky = mxGetPr(mx_ky);
	const double *r_dat = mxGetPr(mx_datr);
	const double *i_dat = mxGetPr(mx_dati);
	const double *p_kernel = mxGetPr(mx_kernel);

	double *r_fm, *i_fm;

	if (ndat1 > ndat2)
		ndat = ndat1;
	else
		ndat = ndat2;

	nn = nim[0];
	W = www[0];
	/* printf("%d %d %d %d %d\n",nn,nim,ndat,ndim,W);
	   for (i = 0; i < ndat;i++)
	   printf("%f %f\n",*(r_dat + i),*(i_dat + i)); */
	/* create a new array and set the output pointer to it */
	if (N != 1)
	{
		mexFail0("N=1 done only");
	}
	else
		plhs[0] = mxCreateDoubleMatrix(nim[0]*nim[0],1,mxCOMPLEX);
	r_fm = mxGetPr(plhs[0]);
	i_fm = mxGetPi(plhs[0]);
	kbgrid(p_kx,p_ky,ndat,r_dat,i_dat,nn,W,p_kernel,nk,r_fm,i_fm);

	return 1;
}


/* The gateway routine. */
void mexFunction(
		int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
	/* check for the proper number of arguments */
	if (nrhs != 7)
		mexFail("7 inputs required: (kx, ky, datr, dati, nim, W, kernel)");
	if (nlhs > 1)
		mexFail("Less than one output arguments.");

	if (!kbgrid_mex(plhs,
				prhs[0], prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], prhs[6]))
		mexFail("kbgrid_mex failed");

	return;
}
