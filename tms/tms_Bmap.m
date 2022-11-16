function tms_Bmap(xres,names, doRecon)
%function tms_Bmap(xres, Pfiles(cell array), doRecon?....)
%
% (c) Luis Hernandez-Garcia and Sangwoo Lee @ University of Michigan
% maintained by hernan@umich.edu
%
%
% example : 
% P = dirnames('P*')
% tms_Bmap(128,P, 1)
%
% assumes that the first half of the P files contains the images 
% collected while the current was on
%
% reconstruct magnitude images and field maps using gsp21b. 
% corregister them to the first volume
% -> write out the rigid body transformation matrices
% calculate the B vector from its projections onto the Z axis
%
% Needs the FSL library, Doug's gsp recon, and the script "flirt_batch"
% in addition to my image manipulation stuff

FSL=1;  % use either fsl or SPM package for registration/reslicing
SPM=0;
if length(names)<6
    error('the minimum number of measurement is 3 data points.');
end;
if mod(length(names),2)~=0
    error('the measurement should be in pairs');
end;

num_datapoints=length(names)/2;
ONfiles = names(1:num_datapoints);
OFFfiles = names(num_datapoints+1:end);

if doRecon
  %!rm -f P*
  %!rm -f *.txt
  !rm -f *.nii.gz
  !rm -f *.img
  !rm -f *.hdr
  !rm -f ref*
  !rm -f *.mat
  !rm -f vol*
  !rm -f rvol*
  !rm -f sl*
  !rm -rf rtmp.mat*

%%% reconstruct the structural images and store them as analyze format in
%%% 'vol*.img' and 'vol*.hdr' and save them in 01 03 05... fashion
%%%% these are images used to calculate the transformation matrices.
  for ii=num_datapoints:-1:1
    disp(sprintf('reconstructing %dth set of structural images...',ii));
    str = sprintf('!gsp21b -n %d -fx -fy -m %s',xres, char(OFFfiles(ii)))
    eval(str)
    str = sprintf('!gsp21b -A -n %d -fx -fy -h %s',xres, char(OFFfiles(ii)))
    eval(str)
    P = dir('vol_*.img');
    root = P(1).name; oldname = root(1:end-4);
    str = sprintf('!mvimg %s vol_open%02d', oldname, ii*2-1)
    eval(str)
  end;
%%% reconstruct the open and on images and get the field maps and store
%%% them into an analyze format as 'open*.img' and 'open*.hdr'
  for ii=1:num_datapoints
    disp(sprintf('reconstructing %dth set of phase images...',ii));
    mapdel=2500e-6;
    eval(['!gsp21b -0 -p -fx -fy -n ' sprintf('%d ',xres) char(OFFfiles(ii))]); 
    P = dir2names('vol*.hdr');
    aa=char(P(1));
    hdr1=read_hdr(aa);
    hdr1.datatype=16;  % fake the header file
    nsl = hdr1.zdim;
    map=fmap(nsl,xres,mapdel);
    write_img(sprintf('BZ_open%d.img',ii),map,hdr1);
    write_hdr(sprintf('BZ_open%d.hdr',ii),hdr1);


    eval(['!gsp21b -0 -p -fx -fy -n ' sprintf('%d ',xres) char(ONfiles(ii))]); 
    map=fmap(nsl,xres,mapdel);
    write_img(sprintf('BZ_on%d.img',ii),map,hdr1);
    write_hdr(sprintf('BZ_on%d.hdr',ii),hdr1);
  end;
end
% clean up residual files
!rm -f ref*
!rm -f sl*

% realign images to calculate the transformation matrices.  using FSL
%! avwmerge -t tmp vol_open*.img
%! mcflirt -in tmp -out rtmp -refvol 0 -cost normcorr -verbose -stats -plots -mats

%! avwmaths tmp.hdr -thr 800 ttmp.hdr
%! mcflirt -in tmp -out rtmp -refvol 2 -smooth 5 -cost mutualinfo -verbose -stats -plots -mats 
%! avwsplit rtmp.img
%! rm -f rtmp.img rtmp.hdr tmp.hdr tmp.img
! flirt_batch vol_open01 vol_open    

%keyboard; % check registration
smoother('BZ_on',4);
smoother('BZ_open',4);

% applying the transformations to the phase maps
disp('Apply transformations to the phase images ...');
Pnames = dir('s_BZ_open*.img');
Pnames2 = dir('s_BZ_on*.img');
%MATnames = dir('rtmp.mat/MAT*');
MATnames = dir('./MATS/MAT*')


for f=1:length(Pnames)
   str = sprintf('!flirt -in %s -ref %s -applyxfm -init MATS/%s -out r%s',...
   Pnames(f).name, Pnames(f).name, MATnames(f).name, Pnames(f).name)
   eval(str)
   str = sprintf('!flirt -in %s -ref %s -applyxfm -init MATS/%s -out r%s',...
            Pnames2(f).name, Pnames2(f).name, MATnames(f).name, Pnames2(f).name)
   eval(str)
end;

%% after running SPM2, (spm fmri) (or FSL view)
%% read in the realigned field maps for on and open cases
rs_open = read_img_series('rs_BZ_open');
rs_on = read_img_series('rs_BZ_on');
num_pix = size(rs_open,2); 

% load in realigned-resliced magnitude images for mask
names = dir('rvol_open0*.img');
maskim=read_img2(names(1).name);

% generate 3D mask
mask=maskim>max(maskim(:))*0.2;

hdr1=read_hdr('s_BZ_open1.hdr');
write_img('bmask.img',mask,hdr1);
write_hdr('bmask.hdr',hdr1);
write_img('maskim.img',maskim,hdr1);
write_hdr('maskim.hdr',hdr1);

% calculate the induced Bz field from the phase difference map and store
% them into Bz_1, Bz_2 ...
gamma_bar=4258; % Hz / Gauss
Bz_measured = (rs_on-rs_open)/mapdel/gamma_bar*2*pi;
Bz_measured(find(isnan(Bz_measured)))=0;
%save Bzs Bz_1 Bz_2 Bz_3 mask

% read in the transformation matrices
% we are interested in the third row of eac transformation matrix
cd MATS
mfiles = dir('MAT_*');
M = zeros(length(mfiles),3);
for ii=1:length(mfiles)
    tmp =load(mfiles(ii).name);
    itmp=inv(tmp);
    M(ii,:) = itmp(3,1:3);
end;
cd ..
% finally solve the problem.
% solve the inverse problem for Bx and By
% following the notation of Hernandez et al.
iM=pinv(M);
disp(sprintf('condition number of aa is %f',cond(M)));
Bx=zeros(1, size(Bz_measured,2));
By=zeros(1, size(Bz_measured,2));
Bz=zeros(1, size(Bz_measured,2));
dim = size(Bz_measured);

         temp2=Bz_measured;
         temp=iM*temp2;
         Bx=temp(1,:);
         By=temp(2,:);
         Bz=temp(3,:);

B_mag=sqrt(Bx.^2+By.^2+Bz.^2);

hdr1=read_hdr('s_BZ_open1.hdr');
save tms_wkspace 

Bxm=Bx.*mask(:)';
Bym=By.*mask(:)';
Bzm=Bz.*mask(:)';

write_img('Bx.img',Bxm,hdr1);
write_img('By.img',Bym,hdr1);
write_img('Bz.img',Bzm,hdr1);
write_hdr('Bx.hdr',hdr1);
write_hdr('By.hdr',hdr1);
write_hdr('Bz.hdr',hdr1);
save results_masked_clipped Bxm Bym Bzm mask

%%%%%%%  now let's do some calculations about the noise   %%%%
Bnoises=[];
for sigma= 0.1:0.1:2  % let the noise (std. dev) in the BZ map measurement be sigma
	varBZ = sigma^2 * eye(size(iM,2));
	varBhat = iM * varBZ * iM';
	sigmaBhat = sqrt(diag(varBhat));
	Bnoises = [Bnoises sigmaBhat];
end
plot([0.1:0.1:2],Bnoises)
title('Noise Amplification');
xlabel('\sigma_{Bz} measured')
ylabel('\sigma_{Bxyz} estimate');
legend('\sigma_{Bx}','\sigma_{By}','\sigma_{Bz}')
fatlines
dofontsize(16)


%%%%%%%%%%%%%%% calculate what the movement parameters were:
units = [180/pi 180/pi 180/pi 1 1 1];

m0 = load('MAT_0000'); p0=getMovParms(m0) .* units;
m1 = load('MAT_0001'); p1=getMovParms(m1) .* units;
m2 = load('MAT_0002'); p2=getMovParms(m2) .* units;
m3 = load('MAT_0003'); p3=getMovParms(m3) .* units;

[p0 ; p1 ; p2; p3]
%keyboard;
%orthoq('B','maskim',[-30 30]);
%keyboard;
if 0
    %% after rotating Bx,By,Bz around manually using Bz in spm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !cp Bz.mat Bx.mat
    !cp Bz.mat By.mat
    !cp Bz.mat maskim.mat
    !cp Bz.img Bogus.img
    !cp Bz.hdr Bogus.hdr

    P1 = dir2names('B?.img');
    P2 = dir2names('maskim*.img');
    P1(end+1)=P2;
    spm_reslice(P1)



    rBx=read_img2([],'rBx');
    rBy=read_img2([],'rBy');
    rBz=read_img2([],'rBz');
    rmaskim=read_img2([],'rmaskim');

    %%% run orthoq
    orthoq('rB','rmaskim',[-30 30]);
end;

function Conditions
rot = 0:0.01:pi/4;
MM = zeros(4,3);
conds = zeros(size(rot));

for r = 1:length(rot)
    parms = [rot(r) 0 0 0 0 0];
    M=makeAffine(parms);
    MM(1,:) = M(3,1:3);
    parms = [0 rot(r) 0 0 0 0];
    M=makeAffine(parms);
    MM(2,:) = M(3,1:3);
    parms = [0 0 rot(r) 0 0 0];
    M=makeAffine(parms);
    MM(3,:) = M(3,1:3);
    MM(4,3)=1
    conds(r) = cond(MM);
end
plot(rot*180/pi,conds)
line([0 40],[10 10])
xlabel('Rotation (degrees)')
ylabel('Condition Number ')
title('Invertability of System')
dofontsize(16)
fatlines

return
