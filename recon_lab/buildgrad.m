% $HeadURL: svn+ssh://klitinas@woo.engin.umich.edu/work/repos/matlab/recon/branches/phase_zero/buildgrad.m $
% $Id: buildgrad.m 1266 2014-03-20 18:12:40Z klitinas $

spup = 1;
opfov = 20/spup;
opxres = 64/spup;
gts = 4e-6;
gslew = 180;
gamp = 2.6; 
nl = 2;
[g,k,t,s]=spiralgradlx3(opfov,opxres,gts,gslew,gamp,nl);
