function eucv=finddist2D(pp,width,height,refno,dmloc);
%
%	finds euclidean distances of the pixels from a distance mapels
%
%	cengizhan 9.16.1996
%
%	Mod 9/1/97--> vectorized version
%
lt=size(pp,2);			% every 2D point in a single column
eucv=zeros(lt,2);
exv=zeros(lt,2);
pp=round(pp);			
%
% 	For the points outside the pixel area find the closest point
%	on the space and make an addition 
%
%for ii=1:lt
%	if (pp(1,ii)<1),
%		exv(ii,1)=1-pp(1,ii);					
%		pp(1,ii)=1;
%	elseif (pp(1,ii)>width),
%		exv(ii,1)=pp(1,ii)-width;
%		pp(1,ii)=width;
%	end;
%	if (pp(2,ii)<1),
%		exv(ii,2)=1-pp(2,ii);
%		pp(2,ii)=1;
%	elseif (pp(2,ii)>height),
%		exv(ii,2)=pp(2,ii)-height;
%		pp(2,ii)=height;
%	end;
%end;
	dum=pp(1,:)<1;
	exv(dum,1)=1-pp(1,dum)';				
	pp(1,dum)=ones(1,nnz(dum));
	
	dum=pp(1,:)>width;
	exv(dum,1)=pp(1,dum)'-width;				
	pp(1,dum)=ones(1,nnz(dum));

	dum=pp(2,:)<1;
	exv(dum,2)=1-pp(2,dum)';				
	pp(2,dum)=ones(1,nnz(dum));

	dum=pp(2,:)>height;
	exv(dum,2)=pp(2,dum)'-height;				
	pp(2,dum)=ones(1,nnz(dum));
%
% find the distances from ndist2D.mex
%
eucv=ndist2D(pp',width,height,(8+dmloc),refno);	
%
eucv=eucv+exv;						% add the extra distances for the points outside voxel space
