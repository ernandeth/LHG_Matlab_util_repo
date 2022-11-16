function gsinterp2(width,height,testno,fname);
%	format: gsinterp2(width,height,testno,fname)
%	
%	this function interpolates a grayscale row image
%	according to the parameters in paramxxx (xxis given in testno)
%
%	First a lookupxxx is generated using transformation parameters
%	which are in paramxxx. The later should be generated 
%	by an initial alignment using ALIGN on the reference 
%	and test sets.
%
%	Lookupxxx maps the reference plane back to the test plane,
%	so every point in the reference plane can be found by interpolating
%	around its corresponding point (given by xlook and ylook)
%
%	Cengizhan
%	7.10.1997
%
%	Ps:1. fname should be the full name of teh grayscale raw image
%	Ps:2. origx and origy are original image sizes
%
 
disp('.......Generating lookup table!');

xlook=zeros(height,width);
ylook=zeros(height,width);
[xlook,ylook] = meshgrid(1:width,1:height);
dumref=[reshape(xlook,1,height*width); reshape(ylook,1,height*width)];
eval(['load param' int2str(testno) ';']);
eval(['CC=CC' int2str(testno) ';']);
Cnew=CC+tt(2:3);
if exist('tt'),
	RR=getRz2D(tt(1));	
	invRR=inv(RR);
	newdumref=getnewP2D(invRR,-tt(2:3),dumref,Cnew);
end;
%		
xlook=reshape(newdumref(1,:),height,width);
ylook=reshape(newdumref(2,:),height,width);		
%   
%
disp('.......Loading the image!');
%
% interpolate    
%
%[img,map] = tiffreadmod(fname);       % put a flag here later to load and save tiff images XXXXXXXXXXXXX

fid=fopen(fname,'rb');
img=fread(fid,'uint8');
img=(reshape(img,width,height))';
fclose(fid);

disp('.......Interpolating!');

img2=zeros(width,height);
XX = 1:width;
YY = (1:height)';
[Xi,Yi] = meshgrid(XX,YY);

img2 = interp2(Xi,Yi,img,xlook,ylook);

% save the image
%
% tiffwritemod(img2,map,[fname '_trans'])
disp('.......Saving the new image!');

fid=fopen(['n' fname],'wb');
fwrite(fid,img2','uint8');
fclose(fid);

%
disp(['finished image ' int2str(testno)]);
%
