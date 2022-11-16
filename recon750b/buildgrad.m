% $Id: buildgrad.m 304 2012-10-24 15:47:14Z klitinas $

spup = 1;
opfov = 20/spup;
opxres = 64/spup;
gts = 4e-6;
gslew = 180;
gamp = 2.6; 
nl = 2;
[g,k,t,s]=spiralgradlx3(opfov,opxres,gts,gslew,gamp,nl);