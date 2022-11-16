clc
clear all

global diagn count_iter % Define global variables for storing all important parameters in the cost function
diagn = struct();
count_iter = 0;

% tstart=tic;

n = 64; % Size of m(x,y) For matrix size = 64 x 64
% n = 16; %For matrix size = 16 x 16
dx = 4; % Pixel size for dx
dy = 4; % Pixel size for dy
gambar=4257; % Hz/G

%---------------------------Memory allocation------------------------------
signal = zeros(1,n*n); % s(t) = sum(sum(m(x,y) exp(-i*phase1))) without external field
signal_AC = zeros(1,n*n); % s(t) = sum(sum(m(x,y) exp(-i*phase1) exp(-i*phase2))) with external field
Bfield_new_AC = zeros(n,n,n*n);% External magnetic field AC

kx = zeros(1,n*n); % Memory allocation for k-space trajectory
ky = zeros(1,n*n);


% Use for phase calculations
% phase_voxel = zeros (n,n,n*n);
% phase_voxel_DC = zeros (n,n,n*n);
% phase_voxel_AC = zeros (n,n,n*n);




%--------------------------------------------------------------------------


%Define the original image

%Square phantom
% m1 = zeros(n,n);
% m1(22:44,22:44) = 1;

%Shepp-Logan Phantom
m1=phantom(n,n);
m = m1(:);

% X and Y cordinates for each voxel
x = linspace(-(n/2),(n/2)-1,n)*dx;
y = linspace(-(n/2),(n/2)-1,n)*dy;
[xx_mat yy_mat] = ndgrid(x,y);
xx = xx_mat(:);
yy = yy_mat(:);
time = linspace(0,32e-3,n*n); % n=64;
% time = linspace(0,8e-3,n*n);% n=16;
dt=time(2)-time(1); %Step in time


%--------------------------------------------------------------------------
%Calculation of Bfield
scaling = 1e3; % Increase the amplitude of the current to 1 A from 1mA

%Define the static magnetic field B0
B0x = zeros(n,n);
B0y = zeros(n,n);
B0z = ones(n,n) * 1.5 .* 1e4;

%Load the original field map from the Biot Savart calculator
load 'wire_figure8coil_100pts.mat'

% biot_3D does calculations in Tesla. All calculation in this code are done
% in gauss. 1 T = 10^4 gauss
Bfx = bx(:,:,7) .* 1e4;
Bfy = by(:,:,7) .* 1e4;
Bfz = bz(:,:,7) .* 1e4;

% Bfield = abs(B0 + B2) - abs(B0), where B2= DC field map
babs = sqrt( (Bfx + B0x) .^ 2 + (Bfy + B0y) .^ 2 + (Bfz + B0z) .^ 2 );
B0_mag = sqrt( B0x .^ 2 + B0y .^ 2 + B0z .^ 2 );

Bfield = (babs - B0_mag ) .* scaling;


% Define model for AC magnetic field
freq=1e4; % Frequency
amp = 1e-3; % Manipulate amplitude according to requirements, If amp=1 current = 1A
bf_func = @(t) Bfield .* cos(2*pi*freq*t ) .* amp;

% Store the Bfield evolving as a function of time
Bfield_new_AC(:,:,1) = bf_func(time(1));


%--------------------------------------------------------------------------
% K-space trajectory caluclations

dkx = 1/(n*dx); % Step in k-space along x direction
dky = 1/(n*dy);  % Step in k-space along y direction

FOVx = 1/dkx; % Field of view along x
FOVy = 1/dky; % Field of view along y

Gx = 0.1175;
Gy = 0.1175;

% Starting position in k-space
kx(1) = - 0.25e-3 * Gx * gambar  ;
ky(1) = - 0.25e-3 * Gy * gambar ;


%--------------------------------------------------------------------------


Bf_AC = Bfield_new_AC(:,:,1);
Bf_AC = Bf_AC(:);

%This loop does calculation for each voxel for only the first time point.
%Have to remove this lopp and incorporate it in the main loop
for r = 1 :  n*n
    
    signal_temp = (m(r) * exp(-1i * (2 * pi * (kx(1) * xx(r) + ky(1) * yy(r)))));
    
    signal(1) = signal_temp + signal (1); % Signal without external magnetic field
    
    
    signal_temp_AC= m(r) * exp(-1i * (2 * pi * (kx(1) * xx(r) + ky(1) * yy(r)))) * exp(-1i * gambar * Bf_AC(r) );
    
    signal_AC(1) = signal_temp_AC + signal_AC(1); % Signal with external magnetic field
    
end


count = 1; % counter for changing directiong in k-space after 'n' points

for t = 2: n*n
    fprintf(' Sample %g of %g \n',t,n*n);
    
    count = count + 1; % Flip after 'n' points in k-space
    itime = t; % Counter for time
    Gy=0;
    
    %----------------------------------------------------------------------
    % Calculation for kx
    
    if ( count == (n+1))
        Gx=-Gx;
        kx(itime) = kx(itime-1);
    else
        temp = gambar * Gx * time(itime) / (itime-1);
        kx( itime) =  kx(itime-1) + temp;
    end
    
    %----------------------------------------------------------------------
    % Calculation for ky
    
    
    temp = gambar * Gy * time(itime) / (itime - 1);
    ky( itime) = ky(itime-1) + temp;
    
    if (count == (n+1))
        %                 Gy = 0.0092*4; % for n = 16;
        Gy = 0.0092; % for n = 64;
        count = 1;
        ky(itime ) = ky(itime - 1) + gambar * 0.1e-3 * Gy ;
        
    end
    
    %-----------------------------------Calculation for Bfield-------------
    %  Bfield evolving with time.
    Bfield_new_AC(:,:,itime) = (bf_func(time(itime))) + Bfield_new_AC(:,:,itime-1);
    Bf_AC = Bfield_new_AC(:,:,itime);
    Bf_AC = Bf_AC(:);
    
    
    %----------------------------------------------------------------------
    
    for r = 1 :  n*n
        
        
        %-----------------Without magnetic field--------------------------
        signal_temp = (m(r) * exp(-1i * (2 * pi * (kx(itime) * xx(r) + ky(itime) * yy(r)))));
        
        signal(itime) = signal_temp + signal (itime);
        
        %         phase_temp(r) = phase;
        
        
        %-----------------With AC field------------------------------------
        signal_temp_AC= m(r) * exp(-1i * (2 * pi * (kx(itime) * xx(r) + ky(itime) * yy(r)))) * exp(-1i * gambar * Bf_AC(r) );
        
        signal_AC(itime) = signal_temp_AC + signal_AC(itime);
        
        %         phase_temp_AC(r) = phase + phase_Bfield_AC;
    end
    
    % Phase calculations
    
    %     phase_voxel(:,:,itime) = phase_temp;
    %     phase_voxel_DC(:,:,itime) = phase_temp_DC;
    %     phase_voxel_AC(:,:,itime) = phase_temp_AC;
    %
    
    
end


%Parameters needed to be passed to cost function stored as one structure

parms.m = m;
parms.t = time;
parms.size = n;
parms.xcord = xx;
parms.ycord = yy;
parms.gambar = gambar;
parms.freq = freq;
parms.m1 = m1;
parms.kx = kx;
parms.ky = ky;
parms.time =  time;
parms.amp = amp;


% Reshape the signal onto the k-space trajectory
kspace = reshape(signal,n,n);
kspace(:, 2:2:end) = kspace( end:-1:1, 2:2:end);
kspace = flipud(kspace);

kspace_AC = reshape(signal_AC,n,n);
kspace_AC(:, 2:2:end) = kspace_AC( end:-1:1, 2:2:end);
kspace_AC = flipud(kspace_AC);

image_1 = ifftshift(ifft2(fftshift(kspace)));
image_3 = ifftshift(ifft2(fftshift(kspace_AC)));

% prog_time = toc(tstart)


% Display figures
% figure(1),subplot(131),imagesc(x,y,m1),axis image,colormap(gray),xlabel('x(mm)'),ylabel('y (mm)'),title('Original Image');
% figure(1),subplot(132),imagesc(linspace(-0.2,.2,64),linspace(-0.2,0.2,64),abs(log(1+kspace))),axis image,xlabel('kx (1/mm)'),ylabel('ky (1/mm)'),colormap(gray),title('K-space');
% figure(1),subplot(133),imagesc(y,x,abs(image_1)),axis xy,axis image,xlabel('x(mm)'),ylabel('y (mm)'),colormap(gray),title('Reconstructed image');
%
% figure(2),subplot(131),imagesc(x,y,m1),axis image,colormap(gray),xlabel('x(mm)'),ylabel('y (mm)'),title('Original Image');
% figure(2),subplot(132),imagesc(linspace(-0.2,.2,64),linspace(-0.2,0.2,64),abs(log(1+kspace_AC))),axis image,xlabel('kx (1/mm)'),ylabel('ky (1/mm)'),colormap(gray),title('K-space');
% figure(2),subplot(133),imagesc(x,y,abs(image_3)),axis xy,axis image,xlabel('x(mm)'),ylabel('y (mm)'),colormap(gray),title('Distorted image');


% Phase calculations
% xcord = 34;
% ycord = 41;
% ans(1,:) = phase_voxel_AC(xcord,ycord,:);
% figure;plot(time,ans)

%Starting estimate for lsqnonlin
Bfield_DC_map = zeros(n,n);
est0 = Bfield_DC_map;

optvar=optimset('lsqnonlin');
optvar.Display = 'iter';
optvar.MaxIter = 400;
optvar.MaxFunEvals = 4097*400;

[est res] = lsqnonlin(@signal_func,est0,...
    [],[],...
    optvar,...
    parms,...
    signal_AC);

residual = zeros(1,length(diagn));

for ii = 1:length(diagn)
    residual(ii) = diagn(ii).res;
end

% % % figure;plot(1:length(diagn),residual);

% save 'Experimentresult_phantom.mat' est res Bfield diagn;


