

%% combine individual files into mag and phase sequences; 
% phase images are phase unwrapped using 'pm_unwrap()'
% Folder '170320phantom_tdcs' saves original nii files
% Sub-folder '0320' saves mag and phase
addpath(genpath(pwd));

mag = zeros(14,64*64*60);% magnitude
for n=1:7 % 7 positions
    filename = ['./data/170320phantom_tdcs/vol_P_pos',num2str(n),'_off.nii'];
    [tmp, tmp_h] = read_nii_img(filename);
    tmp = reshape(tmp,1,64*64*60);
    mag(2*n-1,:) = tmp; % 1,3,5... odd numbers are current off
    filename = ['./data/170320phantom_tdcs/vol_P_pos',num2str(n),'_2ma.nii'];
    [tmp, tmp_h] = read_nii_img(filename);
    tmp = reshape(tmp,1,64*64*60);
    mag(2*n,:) = tmp;   % 2,4,6... even numbers are current on
end
head_14 = tmp_h;
head_14.dim(5) = 14;
write_nii('./data/170320phantom_tdcs/0320/mag_14.nii',mag,head_14,0);

phs = zeros(14,64*64*60);% phase
for n=1:7 % 7 positions
    filename = ['./data/170320phantom_tdcs/p_vol_P_pos',num2str(n),'_off.nii'];
    [tmp, tmp_h] = read_nii_img(filename);
    tmp = reshape(tmp,1,64*64*60);
    phs(2*n-1,:) = tmp;  % 1,3,5... odd numbers are current off
    filename = ['./data/170320phantom_tdcs/p_vol_P_pos',num2str(n),'_2ma.nii'];
    [tmp, tmp_h] = read_nii_img(filename);
    tmp = reshape(tmp,1,64*64*60);
    phs(2*n,:) = tmp;    % 2,4,6... even numbers are current on
end

%unwrap phase
mag(mag<1000)=0; % mask
i=sqrt(-1);
raw = mag .* exp(i*phs/1000); % complex image
V_x = 64;
V_y = 64;
V_z = 60;

for n=1:14
    tmp = reshape(raw(n,:),V_x,V_y,V_z);
    tmp = 1000*pm_unwrap(tmp,[V_x,V_y,V_z]); % phase unwrap
    phs(n,:) = reshape(tmp,1,64*64*60);
end

head_14 = tmp_h;
head_14.dim(5) = 14;
write_nii('./data/170320phantom_tdcs/0320/phs_14.nii',phs,head_14,0);

%% positions realignment
% first use 'spm_coreg()' to obtain rotate matrix (current on files only, reduce running time)
% then realign mag and phase
% the rotate matrix is saved into 'mat' file, for invert matrix

addpath(genpath(pwd));
V_m = spm_vol(fullfile(cd,'./data/170320phantom_tdcs/0320/mag_14.nii'));
V_p = spm_vol(fullfile(cd,'./data/170320phantom_tdcs/0320/phs_14.nii'));

Vref = V_m(1,1); % set first position as reference
Vtgt = V_m;
% realign
for n=2:2:14 % using current on images to find rotate matrix
    x = spm_coreg(Vref, Vtgt(n)); % find rotate matrix
    mat_data(n).mat = spm_matrix(x);
end
save mat_data_ref_phantom_0320 mat_data

load('mat_data_ref_phantom_0320.mat')
% reslice mag
Vtgt = V_m;
for n=2:2:14
    mat = mat_data(n).mat;
    xform_m = inv(mat)*Vtgt(n).mat ;
    Vtgt(n-1).mat=xform_m; % current on and off are same position, share same rotate matrix
    spm_reslice([Vref, Vtgt(n-1)],struct('which',1));%to change data, first need Vtgt.mat=xform;
    Vtgt(n).mat=xform_m;
    spm_reslice([Vref, Vtgt(n)],struct('which',1));%to change data, first need Vtgt.mat=xform;
end

% reslice phase
Vtgt = V_p;
for n=2:2:14
    mat = mat_data(n).mat;
    xform_m = inv(mat)*Vtgt(n).mat ;
    Vtgt(n-1).mat=xform_m; % current on and off are same position, share same rotate matrix
    spm_reslice([Vref, Vtgt(n-1)],struct('which',1));%to change data, first need Vtgt.mat=xform;
    Vtgt(n).mat=xform_m;
    spm_reslice([Vref, Vtgt(n)],struct('which',1));%to change data, first need Vtgt.mat=xform;
end

%% combine complex images; divide phase difference; obtain bz and b_xyz
% 'raw' is the complex data combined from mag and phase
% phase difference is converted into magnetic => bz
% b_xyz is the inverted original magnetic with x y z projections, inverted from bz
addpath(genpath(pwd));
load('mat_data_ref_phantom_0320.mat')

%complex trans
[mag, mag_h] = read_nii_img('./data/170320phantom_tdcs/0320/rmag_14.nii');% read_nii_img gives both data and head
[phs, phs_h] = read_nii_img('./data/170320phantom_tdcs/0320/rphs_14.nii');
i=sqrt(-1);
%mask
mag(mag<3000)=0;
raw = mag .* exp(i*phs/1000); % complex images

%3D smoothing
for n=1:size(mag,1)
    tmp = reshape(raw(n,:),64, 64, 60);
    %tmp = smooth3(tmp,'gaussian',7); % ordinary 3D smooth
    tmp = mrfilter(tmp,[2.8125,2.8125,4]); % auto kernel size 3D smooth
    tmp = reshape(tmp,1,64*64*60);
    raw(n,:) = tmp;
end

%save real and imaginary part of complex images
p_raw=phs_h;
p_raw.dim(5)=14;
write_nii('./data/170320phantom_tdcs/0320/nii_raw_real.nii',real(raw),p_raw,0);
write_nii('./data/170320phantom_tdcs/0320/nii_raw_imag.nii',imag(raw),p_raw,0);

TE = 5e-3;%35 ms
GAMMA = 42.56 * 1e6 / 1e4;   % gyromagentic constant in Hz/Gauss
F = 1/(GAMMA*2*pi*TE); 

%divide phase difference, multiply scalar to get magnetic
bz = zeros(7,length(raw));
for n=1:7
    tmp = F*angle(raw(2*n,:)./raw(2*n-1,:));
    tmp(isnan(tmp))=0;
    bz(n,:) = tmp;
end

%rotation transform matrix
u = zeros(size(bz,1),3);
for n=1:length(u)
    u(n,:)=mat_data(2*n).mat(3,1:3);
end
u_inv = pinv(u); %inverse rotation transform matrix

%b_xyz: inverted original bx by bz
b_xyz = u_inv*bz;
p_b_xyz=phs_h; % head settings
p_b_xyz.dim(5)=3;
write_nii('./data/170320phantom_tdcs/0320/nii_b_xyz_unwrap.nii',b_xyz*10^3,p_b_xyz,0);
%nano Gauss

p_bz=phs_h;
p_bz.dim(5)=length(u);
write_nii('./data/170320phantom_tdcs/0320/nii_bz_unwrap.nii',bz*10^3,p_bz,0);

%% interpolation to get cubic voxel; current calculation
% bz_q is the interpolated bz, unit of micro Gauss.
% b_xyz here is inverted from bz_q

x = [1:64]*2.8125;
y = [1:64]*2.8125;
z = [1:60]*4;
[X,Y,Z] = meshgrid(x,y,z);

x = [1:64]*2.8125;
y = [1:64]*2.8125;
z = [1:86]*2.8125;
[Xq,Yq,Zq] = meshgrid(x,y,z);

bz_q = zeros(7,64*64*86);
for n=1:7
    tmp = reshape(bz(n,:),64, 64, 60) *10^6;%multiply 10^6, micro Gauss
    tmp = interp3(X,Y,Z,tmp,Xq,Yq,Zq); % interpolate from 64*64*60 to 64*64*86
    tmp(isnan(tmp))=0;
    bz_q(n,:) = reshape(tmp,1,64*64*86);
end

p_bz=phs_h; % head settings
p_bz.dim(5)=length(u);
p_bz.dim(4)=length(z);
p_bz.pixdim(4)=2.8125;

write_nii('./data/170320phantom_tdcs/0320/nii_bz_unwrap_interp.nii',bz_q,p_bz,0);%micro Gauss

%interpolated b_xyz
b_xyz = u_inv*bz_q;
p_b_xyz=phs_h; % head settings
p_b_xyz.dim(5)=3;
p_b_xyz.dim(4)=length(z);
p_b_xyz.pixdim(4)=2.8125; % unit=mm
write_nii('./data/170320phantom_tdcs/0320/nii_b_xyz_unwrap_interp.nii',b_xyz,p_b_xyz,0);%micro Gauss

%current calculation
bx = reshape(b_xyz(1,:),64, 64, 86);
by = reshape(b_xyz(2,:),64, 64, 86);
bz = reshape(b_xyz(3,:),64, 64, 86);
tmp=[];
%[curlx,curly,curlz]=curl(Xq,Yq,Zq,bx,by,bz);
u0 = 4*pi*10^(-7); % unit=H/m. Permeability (electromagnetism)
dx = 2.8125*10^(-3); % unit=m
scalar = 10^(-4) /(u0*dx) *10^(-4); % first 10^(-4) is for Gauss to T, second 10^(-4) fit nii file value range

[curlx,curly,curlz] = curl(bx,by,bz); % curl function

tmp(1,:) = reshape(curlx,1,64*64*86) *scalar; %unit=uA/mm2*(10^-2), every '1' in nii is 10^-8*(A/mm2)
tmp(2,:) = reshape(curly,1,64*64*86) *scalar; %using such scalar to fit nii file value range
tmp(3,:) = reshape(curlz,1,64*64*86) *scalar; 
write_nii('./data/170320phantom_tdcs/0320/nii_current.nii',tmp,p_b_xyz,0);

tmp2 = tmp(1,:).^2+tmp(2,:).^2+tmp(3,:).^2;
tmp2 = tmp2.^0.5; % magnitude of current
p_b_xyz.dim(5)=1;
write_nii('./data/170320phantom_tdcs/0320/nii_current_mag_test.nii',tmp2,p_b_xyz,0);





