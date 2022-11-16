function m3Daff(testno,refno,width,height,depth,saveparam,randsamp,exitflag,plotflag,plotfinal,remout,afftresh,maxiter,alfa,vv,dmloc,maxiterloop);
%
%	m3Daff(testno,refno,width,height,depth,saveparam,randsamp,...
%		exitflag,plotflag,plotfinal,remout,afftresh,maxiter,alfa,vv,dmloc)
%	returns no official output;
%		- testno is an integer (eg: if it is 1 this program loads test1)
%		- refno is an integer (eg: if it is 1 this program loads ref1)
%		- width,height,depth are dimensions of voxel space
%		- if saveparam==1 it saves param1.mat(if testno is 1)
%			  which contains the transformation parameters R:
%				(3 by 3 matrix) 
%		  it also save err1.mat (if testno is 1) which contains err:
%			 each row gives mean, std, and max errors at each iteration
%				regardless you have chosen a sample or not
%				last row of err1.mat is the map errors of ALL points
%
%		- if randsamp==0, than all of the points are used, else
%			randsamp points are chosen randomly
%		- if exitflag==0, distancemap is kept at memory(usefull for
%			back to back alignments using the same reference),
%		    if exitflag==1, it is flushed
%		- if plotflag is 1, you can follow a sample of
%			points move as they are optimized (not recommended !)
%		- if plotfinal is 1, you can see the alignment after 
%			they are optimized
%		- if dmloc is 0, DM is loaded and used from RAM
%			if dmloc is 1, DM is read from disk as needed
%
%	It uses  uses ndist9.mex, finddist2.m 
%			ndist9.mex file is macintosh type native(either 68k or PPC)
%			it has to be recompiled from ndist9.c for different machines
%
%	OUTLINE:
%	The program first looks whether distancemapxx (DM) is in memory or not.
%	If it can not find it in RAM it looks inside the current folder 
%	for a file named 'distancemapxx'  (xx is refno)
%	if it finds distancemapxx, it is loaded
%	if there is no distancemapxx, refxx.mat (which contain ref) is loaded
%		and it's distancemapxx is generated AND saved
%	than it loads desired testxx.mat (xx is testno !!)
%		which should contain variable test
%	both test and ref are in 3 column format where each point is in a row
%
%
%	Good Luck !!	
%	
%	Cengizhan Ozturk
%	10.25.96
%	ozturkc@dunx1.ocs.drexel.edu
%	
%  PS: Other parameters for optimization and STOP criteria:
%						
%	afftresh, limit in affine parameter change 
%	maxiter, maximum no of iterations 
%	alfa,	optimization parameter 1
%	vv, optimization parameter 2 (divider) 
%	maxiterloop, maximum case 3 loop no
%
dim=3;							
%=================================================
%   IF YOU CHANGE ANYTHING BELOW THIS LINE
%	MAKE A NOTE AND BASICALLY YOU ARE ON YOUR OWN
%=================================================
%
% Check for distance map, load the distance map if none is present
%
eval(['global dist' int2str(refno) ';']);
eval(['pdist=dist' int2str(refno) ';']);
disp('.......Looking for distance map!');
if isempty(pdist),
	fid=fopen(['distancemap' int2str(refno)] ,'r+');		
	if (fid>0)&(dmloc<1),
		disp('.......Found distance map in the folder!');
		disp('.......Loading distance map!');
		st=ndist9([0 0 0],width,height,depth,6,refno);		
		disp('.......Starting optimization!');
		eval(['dist' int2str(refno) '=1;']);
		fclose(fid);
	end;
	if (fid==-1),
		disp('.......Generating distance map!');
		eval(['load ref' int2str(refno) ';']);
		if dmloc==0,
			st=ndist9(round(ref),width,height,depth,2,refno);    
			disp('.......Saving distance map!');
			st=ndist9([0 0 0],width,height,depth,4,refno);
			eval(['dist' int2str(refno) '=1;']);
		elseif dmloc==1,
			st=ndist9(round(ref),width,height,depth,3,refno);
		end;
		disp('.......Starting Optimization!');
	end;	
else,
	disp('....Found the distancemap in memory, starting optimization!');
end;
%
%
% initial parameters 
%
R=eye(dim); 				
%
% load the test.data																		Add the generation of the data here;
%
eval(['load test' int2str(testno) ';']);
NN=size(test,1);
%
Ctest=mean(test);							% Mod v1.1
test(:,1)=test(:,1)-Ctest(1);				% Mod v1.1
test(:,2)=test(:,2)-Ctest(2);				% Mod v1.1
test(:,3)=test(:,3)-Ctest(3);				% Mod v1.1
if ~exist('ref'), 
	eval(['load ref' int2str(refno) ';']); 	% Mod v1.1
end;	
Cref=mean(ref);								% Mod v1.1
ref(:,1)=ref(:,1)-Cref(1);					% Mod v1.1
ref(:,2)=ref(:,2)-Cref(2);					% Mod v1.1	
ref(:,3)=ref(:,3)-Cref(3);					% Mod v1.1	
%
if randsamp,
	N=randsamp;						% choose a sample if desired
	P=zeros(N,dim);					 
	samp=rand(N,1)*(NN-1);
	samp=floor(samp+1);
	P=test(samp(:),1:dim);
else,
	Nt=20;
	samp=rand(Nt,1)*(NN-1);
	samp=floor(samp+1);
	P=zeros(NN,dim);
	N=NN;
	P=test;
end;	
%
%
P=P';									% every point in a column from now on 
%
%  DRAWING ROUTINES
%
if plotflag,
	plotfancy(test);
	hold;
	if ~exist('ref'), 
		eval(['load ref' int2str(refno) ';']);	
	end;
	plot3(ref(:,1),ref(:,2),ref(:,3),'r.');
	if randsamp,
		plot3(P(1,:),P(2,:),P(3,:),'o');
	end;
	drawnow;
end;

NewP=zeros(dim,N);
PrevP=zeros(dim,N);
NewP=P;
PrevP=P;													

v=zeros(N,dim);
d=zeros(N,1);
v=finddist2(NewP,width,height,depth,refno,dmloc);		
d=sqrt((v.*v)*[1;1;1]);	
err(1,1)=mean(d);
err(1,2)=std(d);
err(1,3)=median(d);
err(1,4)=max(d);
							
e=sum(d);
disp(['Initial error is:' num2str(e)]);

NewP1=zeros(dim,N);
NewP2=zeros(dim,N);
NewP3=zeros(dim,N);
FinP=zeros(dim,NN);

v1=zeros(N,dim);
d1=zeros(N,1);
v2=zeros(N,dim);
d2=zeros(N,1);
vx=zeros(N,dim);
dx=zeros(N,1);
vf=zeros(NN,dim);
df=zeros(NN,1);

R1=zeros(dim);
R2=zeros(dim);
Rx=zeros(dim);
Rf=zeros(dim);

fla=1;
iter=1;
ind=1;
rempno=1;
indbef=0;
while fla,
	disp(['iter = ' int2str(iter)]);
	if iter > maxiter, break;fla=0; end;
	%
	%show the corresponding point travels with aline between iterations !!
	%
	if plotflag,
		for i=1:N,
			if i~=indbef
				hh(i)=line([PrevP(1,i) NewP(1,i)],[PrevP(2,i) NewP(2,i)],[PrevP(3,i) NewP(3,i)]);
				%	
				%	If a change of point occured 
				%	we do not draw the line !!!
				%		
			end;
		end;
		set(hh(:),'Color',[0 0 1]);
		drawnow;
	end;
	%
	%	find G
	%
	preG=zeros(N,9);
	G=zeros(9,1);
	A=zeros(9,9);
	%
	%
	for i=1:N,
		if d(i), 
			preG(i,1)=v(i,1)*P(1,i)/d(i);		
			preG(i,2)=v(i,2)*P(1,i)/d(i); 				
			preG(i,3)=v(i,3)*P(1,i)/d(i);
			preG(i,4)=v(i,1)*P(2,i)/d(i);
			preG(i,5)=v(i,2)*P(2,i)/d(i);
			preG(i,6)=v(i,3)*P(2,i)/d(i);
			preG(i,7)=v(i,1)*P(3,i)/d(i);
			preG(i,8)=v(i,2)*P(3,i)/d(i);
			preG(i,9)=v(i,3)*P(3,i)/d(i);
		end;		
	end;
	G(1)=2*sum(preG(:,1).*d(:));
	G(2)=2*sum(preG(:,2).*d(:));
	G(3)=2*sum(preG(:,3).*d(:));
	G(4)=2*sum(preG(:,4).*d(:));
	G(5)=2*sum(preG(:,5).*d(:));
	G(6)=2*sum(preG(:,6).*d(:));
	G(7)=2*sum(preG(:,7).*d(:));
	G(8)=2*sum(preG(:,8).*d(:));
	G(9)=2*sum(preG(:,9).*d(:));
	%
	% find A - a little help from the symmetry !!
	%
	for i=1:9,
		for j=i:9,
			A(i,j)=2*sum(preG(:,i).*preG(:,j));
		end;
	end;
	for i=2:9,
		for j=1:i-1,
			A(i,j)=A(j,i);
		end;
	end;
	%
	PrevP=NewP;
	%
	%	find del1 del2
	%
	del1=inv(A+alfa*eye(9))*(-G);					
	del2=inv(A+(alfa/vv)*eye(9))*(-G);			
	%
	if (all(abs(del1(1:9))<afftresh)), 
		break;fla=0; 
	end;
	%
	R1=R+reshape(del1,dim,dim);
	%R1=R1/sqrt(abs(det(R1)));						%	if det =1
	NewP1=R1*P;
	%
	NewP1x(1,:)=NewP1(1,:)+Cref(1);				% Mod v1.1
	NewP1x(2,:)=NewP1(2,:)+Cref(2);				% Mod v1.1	
	NewP1x(3,:)=NewP1(3,:)+Cref(3);				% Mod v1.1	
	%
	v1=finddist2(NewP1x,width,height,depth,refno,dmloc);	% Mod v1.1	
	d1=sqrt((v1.*v1)*[1;1;1]);
	e1=sum(d1);
	%
	R2=R+reshape(del2,dim,dim);
	%R2=R2/sqrt(abs(det(R2)));  					%	if det =1
	NewP2=R2*P;
	%
	NewP2x(1,:)=NewP2(1,:)+Cref(1);				% Mod v1.1
	NewP2x(2,:)=NewP2(2,:)+Cref(2);				% Mod v1.1	
	NewP2x(3,:)=NewP2(3,:)+Cref(3);				% Mod v1.1	
	%
	v2=finddist2(NewP2x,width,height,depth,refno,dmloc);	% Mod v1.1	
	d2=sqrt((v2.*v2)*[1;1;1]);
	e2=sum(d2);
	%
	%
	if e2<e,
		alfa=alfa/vv;
		R=R2;
		v=v2;
		d=d2;
		e=e2;
		NewP=NewP2;
		disp(['case 1 : ' num2str(e)]);
	elseif (e2>=e)&(e1<e),
		R=R1;
		v=v1;
		d=d1;
		e=e1;
		NewP=NewP1;
		disp(['case 2 : ' num2str(e)]);
	else,
		fla1=1;iter2=1;
		while fla1,
			inalfa=alfa;
			alfa=alfa*vv;
			delx=inv(A+alfa*eye(9))*(-G);				
			Rx=R+reshape(delx,dim,dim);
			%Rx=Rx/sqrt(abs(det(Rx))); %	if det =1
			NewP3=Rx*P;
			%
			NewP3x(1,:)=NewP3(1,:)+Cref(1);				% Mod v1.1
			NewP3x(2,:)=NewP3(2,:)+Cref(2);				% Mod v1.1	
			NewP3x(3,:)=NewP3(3,:)+Cref(3);				% Mod v1.1	
			%
			vx=finddist2(NewP3x,width,height,depth,refno,dmloc);	% Mod v1.1
			dx=sqrt((vx.*vx)*[1;1;1]);
			ex=sum(dx);
			if ex<e, break; fla1=0; end;
			disp([' Looping in case 3 : ' num2str(ex)]);
			iter2=iter2+1;
			if iter2>maxiterloop,fla1=2;alfa =inalfa/vv;break; end;
		end;
		if fla1<2,
			R=Rx;
			v=vx;
			d=dx;
			e=ex;
			NewP=NewP3;
		end;
		disp(['case 3 : ' num2str(e)]);
	end;
	%
	%  preliminary routine to remove the outliars 
	%  during the alignment
	%
	if remout,
		% find the maximum distance point
		indbef=ind;
		[dum ind]=max(d);
		% plot the location 
		if plotflag,
			plot3(NewP(1,ind),NewP(2,ind),NewP(3,ind),'*');
			drawnow;
		end;
		remP(rempno,1)=ind;
		samp2=rand(1)*(NN-1);
		samp2=floor(samp2+1);
		P(1:dim,ind)=test(samp2(:),1:dim)';
		% 
		%v(ind,1:3)=finddist2(test(samp2(:),1:3)',width,height,depth,refno,dmloc);
		%d(ind)=sqrt(v(ind,1:3)*v(ind,1:3)');
		%e=sum(d);
		rempno=rempno+1;
	end;
	iter=iter+1;
	err(iter,1)=mean(d);
	err(iter,2)=std(d);
	err(iter,3)=median(d);
	err(iter,4)=max(d);
end;
%
%
FinP=R*(test');
vf=finddist2(FinP,width,height,depth,refno,dmloc);
df=sqrt((vf.*vf)*[1;1;1]);
iter=iter+1;
err(iter,1)=mean(df);
err(iter,2)=std(df);
err(iter,3)=median(df);
err(iter,4)=max(df);
if plotflag,
	plot3(FinP(1,:),FinP(2,:),FinP(3,:),'b.');
	plot3(NewP(1,:),NewP(2,:),NewP(3,:),'o');
end;
if saveparam, 
	eval(['save param' int2str(testno) ' R Ctest Cref;']);
	eval(['save error' int2str(testno) ' err;']);
end;
%
if plotfinal,
	plotfancy(test);hold;
	if ~exist('ref'), 
		eval(['load ref' int2str(refno) ';']);		
	end;
	plot3(ref(:,1),ref(:,2),ref(:,3),'r.');
	plot3(FinP(1,:),FinP(2,:),FinP(3,:),'g.');
	drawnow;
	%plot3(test(samp(:),1),test(samp(:),2),test(samp(:),3),'o');
	%plot3(FinP(1,samp(:)),FinP(2,samp(:)),FinP(3,samp(:)),'o');
end;
%
disp(['...No of sample points  ' num2str(N)]);
if remout,
	disp(['...No of replaced points during iteration ' num2str(rempno-1)]);
end;
disp(['...Mean map error is  ' num2str(mean(d))]);
disp(['...Std of map error is  ' num2str(std(d))]);
disp(['...Median map error is  ' num2str(median(d))]);
disp(['...Max of map error is  ' num2str(max(d))]);
if randsamp,
	disp(['...Mean map error of all points is ' num2str(mean(df))]);
	disp(['...Std of map error of all points is   ' num2str(std(df))]);
	disp(['...Median map error of all points is  ' num2str(median(df))]);
	disp(['...Max. map error of all points is  ' num2str(max(df))]);
end;
disp('...Calculated affine parameters:'); R

if (exitflag)&(dmloc<1), 
	v=ndist9([0,0,0],width,height,depth,10,refno);
	eval(['dist' int2str(refno) '=[];']);
end;
