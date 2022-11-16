function m2D(testno,refno,width,height,saveparam,randsamp,exitflag,plotflag,plotfinal,remout,angtresh,transtresh,maxiter,alfa,vv,dmloc,maxiterloop);
%
%	m2D(testno,refno,width,height,saveparam,randsamp,exitflag,plotflag,...
%			plotfinal,remout,angtresh,transtresh,maxiter,alfa,vv,dmloc,maxiterloop)
%	returns no official output;
%		- testno is an integer (eg: if it is 1 this program loads test1)
%		- refno is an integer (eg: if it is 1 this program loads ref1)
%		- width,height are dimensions of pixel plane
%		- if saveparam==1 it saves param1.mat(if testno is 1)
%			  which contains the transformation parameters tt:
%				(rot,xtrans,ytrans) 
%			NOTE: rotational parameters are in radian
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
%			points move as they are optimized
%		- if plotfinal is 1, you can see the alignment after 
%			they are optimized
%		- if dmloc is 0, DM is loaded and used from RAM
%			if dmloc is 1, DM is read from disk as needed
%
%	It uses  uses ndist2D.mex, finddist2D.m, getMz2D.m
%			ndist2D.mex file is macintosh type native(either 68k or PPC)
%			it has to be recompiled from ndist2D.c for different machines
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
%	angtresh, limit in rotation								
%	transtresh, limit in translation 
%	maxiter, maximum no of iterations 
%	alfa,	optimization parameter 1
%	vv, optimization parameter 2 (divider)
%	maxiterloop, maximum case 3 loop no  
%
dim=2;pno=3;
%=================================================
%   IF YOU CHANGE ANYTHING BELOW THIS LINE
%	MAKE A NOTE AND BASICALLY YOU ARE ON YOUR OWN
%=================================================
% mod 1: 7.24.1997
%   1. started saving centroid as Cxxx (xxx=testno) in param.mat
%
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
		st=ndist2D([0 0],width,height,6,refno);	
		disp('.......Starting optimization!');
		eval(['dist' int2str(refno) '=1;']);
		fclose(fid);
	end;
	if (fid==-1),
		disp('.......Generating distance map!');
		eval(['load ref' int2str(refno) ';']);
		if dmloc,dmloc=0;end;
		st=ndist2D(round(ref),width,height,2,refno);    
		disp('.......Saving distance map!');
		st=ndist2D([0 0],width,height,4,refno);
		eval(['dist' int2str(refno) '=1;']);		
		disp('.......Starting Optimization!');
	end;
else,
	disp('....Found the distancemap in memory, starting optimization!');
end;
%
% initial parameters zero in this case
% first one is rotation, next two x,y translations
%
tt=zeros(pno,1); 				
%
% load the test.data
%
eval(['load test' int2str(testno) ';']);
NN=size(test,1);
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
% find the centroid of sampled points
%
%  C=mean(P);       MOD MOD MOD v1.1
C=mean(test(:,1:2));
%
P=P';									% every point in a column from now on 
C=C';
%
%  DRAWING ROUTINES
%
if plotflag,
	if ~exist('ref'), 
		eval(['load ref' int2str(refno) ';']);		
	end;
	figure;
	plot(ref(:,1),ref(:,2),'r.');
	hold;
	plot(test(:,1),test(:,2),'g.');
	if randsamp,
		plot(P(1,:),P(2,:),'o');
	end;
	drawnow;
end;
%
% make centroid subtraction	early to save time
%			
dumP=zeros(dim,N);
dumP(1,:)=P(1,:)-C(1);
dumP(2,:)=P(2,:)-C(2);

NewP=zeros(dim,N);
PrevP=zeros(dim,N);
NewP=P;
PrevP=P;													

v=zeros(N,dim);
d=zeros(N,1);
v=finddist2D(NewP,width,height,refno,dmloc);		%	v is vector distances, d is scalar distance
d=sqrt((v.*v)*[1;1]);	
err(1,1)=mean(d);
err(1,2)=std(d);
err(1,3)=median(d);
err(1,4)=max(d);
							
e=sum(d);
disp(['Initial error is : ' num2str(e)]);

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
				hh(i)=line([PrevP(1,i) NewP(1,i)],[PrevP(2,i) NewP(2,i)]);
				%	
				%	A change of point occured therefore
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
	preG=zeros(N,3);
	G=zeros(3,1);
	A=zeros(3,3);
	%
	Mz=getMz2D(tt(1));
	for i=1:N,
		if d(i), 
			preG(i,1)=(v(i,1:2)*(Mz'*dumP(1:2,i)))./d(i);		
			preG(i,2)=v(i,1)./d(i); 				
			preG(i,3)=v(i,2)./d(i);
		end;		
	end;
	G(1)=2*sum(preG(:,1).*d(:));
	G(2)=2*sum(v(:,1));
	G(3)=2*sum(v(:,2));
	%
	% find A - a little help from the symmetry !!
	%
	for i=1:3,
		for j=i:3,
			A(i,j)=2*sum(preG(:,i).*preG(:,j));
		end;
	end;
	for i=2:3,
		for j=1:i-1,
			A(i,j)=A(j,i);
		end;
	end;
	%
	PrevP=NewP;
	%
	%	find del1 del2
	%
	del1=inv(A+alfa*eye(3))*(-G);					
	del2=inv(A+(alfa/vv)*eye(3))*(-G);			
	%
	if (abs(del1(1))<angtresh)&(all(abs(del1(2:3))<transtresh)), 
		break;fla=0; 
	end;
	%
	tt1=tt+del1;
	R1=getRz2D(tt1(1));
	NewP1=getnewP2D(R1,tt1(2:3),P,C);
	v1=finddist2D(NewP1,width,height,refno,dmloc);
	d1=sqrt((v1.*v1)*[1;1]);
	e1=sum(d1);
	%
	tt2=tt+del2;
	R2=getRz2D(tt2(1));
	NewP2=getnewP2D(R2,tt2(2:3),P,C);
	v2=finddist2D(NewP2,width,height,refno,dmloc);
	d2=sqrt((v2.*v2)*[1;1]);
	e2=sum(d2);
	%
	if e2<e,
		alfa=alfa/vv;
		tt=tt2;
		v=v2;
		d=d2;
		e=e2;
		NewP=NewP2;
		disp(['case 1 : ' int2str(e)]);
	elseif (e2>=e)&(e1<e),
		tt=tt1;
		v=v1;
		d=d1;
		e=e1;
		NewP=NewP1;
		disp(['case 2 : ' int2str(e)]);
	else,
		fla1=1;iter2=1;
		while fla1,
			inalfa=alfa;
			alfa=alfa*vv;
			delx=inv(A+alfa*eye(3))*(-G);				
			ttx=tt+delx;
			Rx=getRz2D(ttx(1));
			NewP3=getnewP2D(Rx,ttx(2:3),P,C);
			vx=finddist2D(NewP3,width,height,refno,dmloc);
			dx=sqrt((vx.*vx)*[1;1]);
			ex=sum(dx);
			if ex<e, break; fla1=0; end;
			disp([' Looping in case 3 : ' int2str(ex)]);
			iter2=iter2+1;
			if iter2>maxiterloop,fla1=2;alfa =inalfa/vv;break; end;
		end;
		if fla1<2,
			tt=ttx;
			v=vx;
			d=dx;
			e=ex;
			NewP=NewP3;
		end;
		disp(['case 3 : ' int2str(e)]);
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
			plot(NewP(1,ind),NewP(2,ind),'*');
		end;
		remP(rempno,1)=ind;
		samp2=rand(1)*(NN-1);
		samp2=floor(samp2+1);
		P(1:2,ind)=test(samp2(:),1:2)';
		% 
		%v(ind,1:2)=finddist2D(test(samp2(:),1:2)',width,height,refno,dmloc);
		%d(ind)=sqrt(v(ind,1:2)*v(ind,1:2)');
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
Rf=getRz2D(tt(1));
FinP=getnewP2D(Rf,tt(2:3),test(:,1:2)',C);
vf=finddist2D(FinP,width,height,refno,dmloc);
df=sqrt((vf.*vf)*[1;1]);
iter=iter+1;
err(iter,1)=mean(df);
err(iter,2)=std(df);
err(iter,3)=median(df);
err(iter,4)=max(df);
if plotflag,
	plot(FinP(1,:),FinP(2,:),'b.');
	plot(NewP(1,:),NewP(2,:),'o');
	drawnow;
end;
if saveparam, 
	eval(['CC' int2str(testno) '=C;']);
	eval(['save param' int2str(testno) ' tt CC' int2str(testno) ';']);
	eval(['save error' int2str(testno) ' err;']);
end;
%
if plotfinal,
	figure;
	if ~exist('ref'), 
		eval(['load ref' int2str(testno) ';']);		
	end;
	plot(ref(:,1),ref(:,2),'r.');
	hold;
	plot(test(:,1),test(:,2),'g.');
	plot(FinP(1,:),FinP(2,:),'b.');
	plot(test(samp(:),1),test(samp(:),2),'o');
	plot(FinP(1,samp(:)),FinP(2,samp(:)),'o');
	drawnow;
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
disp(['...Calculated Translation parameters:  ' num2str(-tt(2)) '  ' num2str(-tt(3))]);
disp(['...Calculated Rotation parameters:  ' num2str(-tt(1)/pi*180) ]);

if (exitflag)&(dmloc<1), 
	v=ndist2D([0,0],width,height,10,refno);
	eval(['dist' int2str(refno) '=[];']);
end;
