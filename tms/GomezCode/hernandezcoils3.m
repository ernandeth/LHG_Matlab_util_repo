clear all;

%%%% required information
sourcefolder='/Users/hernan/matlab/tms/CerriSim/cerriV2/';
% conductivitymap='/Users/hernan/matlab/tms/CerriSim/cerriV2/sigma.mat';
% boxndim=[64;64;64];
% boxddim=[.00375;.00375;.005];
conductivitymap='/Users/hernan/matlab/tms/GomezCode/cerriV2/sigmaMap.mat';
boxndim=[256; 256; 124];
boxddim=[.00102;.00102;.00120];

p=6; %%how much zero-padding. p must be >4 and even
flag=1; %%set to 0 to debug 1 to run simulation


% %%%%%%%%%%%Loading brain (independent of code)
load(conductivitymap)
sigma = sigmaMap(1:2:end, 1:2:end, 1:2:end);
boxndim=boxndim/2; boxddim=boxddim*2;
sigmap=sigma;
brainlx=size(sigmap);
sigmap(1:p+brainlx(1),1:p+brainlx(2),1:p+brainlx(3))=0;
sigmap(p/2+1:p/2+brainlx(1),p/2+1:p/2+brainlx(2),p/2+1:p/2+brainlx(3))=sigma;
clear sigma sigmaMap;

boxndim=size(sigmap);
global nx ny nz dx dy dz;

nx=boxndim(1);
ny=boxndim(2);
nz=boxndim(3);

dx=boxddim(1);
dy=boxddim(2);
dz=boxddim(3);

ov02([], sigmap, floor(boxndim(1)/2),floor(boxndim(2)/2),floor(boxndim(3)/2), 0,0,2);
% %%%%%%%%% end loading brain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Coil specs
%%%%%%%%%%%                       Begin Coil specifications
%%%% 1.)The origin of the coordinate system is with respect to the center of the brain
%%%% 2.) an array of coils is defined by 2 two-dimensional cell arrays
%%%% Jvec, rs with indexing as follows:
%%%% rs{coil_index}=(element_id,direction)
%%%% Jvec{coil_index}(element_id,direction)
%%%% coil_index- addresses to a specific coil
%%%% element_id- addresses to a point used to interpolate the coil
%%%% direction= is 1,2,or 3 indicating x,y,or z

% the eta component of the location of the ith interpolant of the jth coil would be
%(rs{j}(i,eta)+rs{j+1}(i,eta))/2

% the eta component current flowing out of the ith interpolant of the jth coil would be
%Jvec{j}(i,eta)

%%% note we do not assume closed currents so in general Jvec has one less
%%% component than rs Thus is it is the responsibility of the user to close
%%% all loops

%%%% element_id- an index for one o
r = 0.03;  % this radius is a fraction of the FOV

aperture = 0;% pi/6;
zpos = 0.5*(nz +p/2)*dz;
R = 0.10;
sfactor = -0.33;
shield_offset = .0625;
theta=0;

current =1e3/1e-4; % this is actually dI/dt in Amps/sec.

%%%%%%%%%%%%%%%%%fCreate coils of choise
[rs{1},Jvec{1}] = make_fig8(current, r, zpos, -theta);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Do not modify This will run the code
hold on;
for i=1:length(rs)
    rs{i}(:,1)=rs{i}(:,1)+nx*dx/2;
    rs{i}(:,2)=rs{i}(:,2)+ny*dy/2;
    rs{i}(:,3)=rs{i}(:,3)+nz*dz/2;
end
if flag==1
    Efield=oneshotcerry(boxndim,boxddim,sigmap,rs,Jvec);
    [Ex1 Ey1 Ez1]=organizefield(Efield);
    save('Efield1_fig8.mat','Ex1','Ey1','Ez1')
else
    seesetup(rs,Jvec,sigmap,ones(size(sigmap)));
end

Eabs = sqrt(Ex1.^2 + Ey1.^2 + Ez1.^2);

for z=1:nz
ov02([], 20*log10(Eabs/max(Eabs(:))), 68, 68,z, 0, -50,1);
pause(0.2);
end



