function eucv=finddist2(pp,width,height,depth,refno,dmloc);
%
%	finds euclidean distances of the voxels
%	within a distancemap
%	uses ndist9.mex
%	cengizhan 11.20.1996
%
%	Mod 9/1/97--> vectorized version
%
lt=size(pp,2);		% every 3D point in a single column in pp
eucv=zeros(lt,3);
exv=zeros(lt,3);
pp=round(pp);			
%
% 	For the points outside the voxel space find the closest point
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
%	if (pp(3,ii)<1),
%		exv(ii,3)=1-pp(3,ii);
%		pp(3,ii)=1;
%	elseif (pp(3,ii)>depth),
%		exv(ii,3)=pp(3,ii)-depth;
%		pp(3,ii)=depth;
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

	dum=pp(3,:)<1;
	exv(dum,3)=1-pp(3,dum)';				
	pp(3,dum)=ones(1,nnz(dum));

	dum=pp(3,:)>depth;
	exv(dum,3)=pp(3,dum)'-depth;				
	pp(3,dum)=ones(1,nnz(dum));

%
% find the distances from ndist9.mex
%
eucv=ndist9(pp',width,height,depth,(8+dmloc),refno);	
%
eucv=eucv+exv;			% add the extra distances for the points outside voxel space

