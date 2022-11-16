function result=k2image(uu,vv,k,ww,N,W)
% function result=k2image(uu,vv,k,ww,N,W)
%
% this function gets gridded k data from non cartesian k data
% then returns the NxN ifft image of it. 
% uu,vv is coordinates of the noncartesian k trajectory data
% ww is the data density precompensation weighting
%	it is obtained by weight_vor, with some modifications on 
%	its behavior at the edges of k-space.

if (nargin<6)
   W=6; % kernel width
end


%disp('new');
xi=-N/2:0.5:N/2-0.5;         % doubles the matrix size and will be
%yi=xi;
yi=-N/2:0.5:N/2-0.5;                       % subsampled after gridding


%Initial range check
if ((max(uu)>N/2) | (max(vv)>N/2) | (min(uu)<(-N/2)) | (min(vv)<(-N/2)))
     sprintf('k-space trajectory is out of range, functionality for k-space deviations outside of -N/2 to N/2 has not been taken into account')
end
    
% fmygrid2 : uses loc and output driven convolution
%loc=getloc(k,N,5,5);
% [temp lx ly bx by]=fmygrid2(uu,vv,k,xi,yi,2.5,2.5,'KB',loc,N);

% mygrid2 : uses input driven convolution
% [temp lx ly bx by]=mygrid2(uu,vv,k,xi,yi,2.5,2.5,'KB');

% ffmygrid : uses output driven convolution
% [temp lx ly bx by]=ffmygrid(uu,vv,k,xi,yi,2.5,2.5,'KB');

% fmygrid3 : uses loc and output driven convolution with weight correction
%loc=getloc(k,N,5,5);
%[temp lx ly bx by]=fmygrid3(uu,vv,k,xi,yi,2.5,2.5,'KB',loc,N);

% ffmygrid4 : uses output driven convolution with weight correction

[temp lx ly bx by]=ffmygrid(uu,vv,ww.*k,xi,yi,W,W,'KB');
% Post gridding weight here.
post_wt = real(ffmygrid(uu,vv,ww.*ones(length(k),1),xi,yi,W,W,'KB'));
l = find(~(post_wt==0));
temp(l) = temp(l)./post_wt(l);

l = find(isnan(temp));
temp(l)=0;

           % to be deleted afterwhile
x=linspace(-0.5,0.5,N*2+1);
x=x(1:2*N);

y=linspace(-0.5,0.5,N*2+1);
y=y(1:2*N);

lx = 2*lx;   %THESE LINES HERE??
ly = 2*ly;

cx = sin(sqrt(pi^2 * lx^2 * x.^2 - bx^2))./ sqrt(pi^2 * lx^2 * x.^2 - bx^2);
cy = sin(sqrt(pi^2 * ly^2 * y.^2 - by^2))./ sqrt(pi^2 * ly^2 * y.^2 - by^2);
[ccx ccy] = meshgrid(cx,cy);
c=ccy.*ccx;
c=c/max(max(c));


result=fftshift(ifft2(fftshift(temp)))./c;


	% to be deleted afterwhile
result=result([N/2+1:3/2*N],[N/2+1:3/2*N]);







