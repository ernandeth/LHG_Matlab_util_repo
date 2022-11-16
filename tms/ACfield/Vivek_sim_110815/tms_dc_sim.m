clear all

% Performs signal equation over space
%%
global diagn count_iter store % Define global variables
diagn = struct();
count_iter = 0;
store = 1;
% 
% filename_fid = '/net/charlie/home/bhatiavr/data/Vnmr_data/data108/gemstrigger_axial_notrigger01.fid/fid';
% filename_procpar = '/net/charlie/home/bhatiavr/data/Vnmr_data/data108/gemstrigger_axial_notrigger01.fid/procpar';


% filename_fid = 'C:\Users\Public\Documents\Data\Vnmr_data\data108\gemstrigger_axial_notrigger01.fid\fid';
% filename_procpar = 'C:\Users\Public\Documents\Data\Vnmr_data\data108\gemstrigger_axial_notrigger01.fid\procpar';
% 
% procpar=struct();
% 
% 
% procpar.gambar = 4257 ; % Gyromagnetic ratio for protons in Hz
% procpar.tr = ReadProcpar('tr',filename_procpar); % Repetition time in seconds
% procpar.te = ReadProcpar('te',filename_procpar); % Echo time in seconds
% procpar.lro = ReadProcpar('lro',filename_procpar); % FOV along readout in cm
% procpar.lpe = ReadProcpar('lpe',filename_procpar); % FOV along pe in cm
% procpar.sw =  ReadProcpar('sw',filename_procpar) ; % Spectral width in Hz
% procpar.at =  ReadProcpar('at',filename_procpar); %Readout time in seconds
% procpar.thk =  ReadProcpar('thk',filename_procpar) * 0.1 ; % Slice thickness in cm
% procpar.gro =  ReadProcpar('gro',filename_procpar) ; % Readout gradient in G/cm
% procpar.gpe =  ReadProcpar('gpe',filename_procpar) ; % Phase Encode gradient in G/cm
% procpar.gror =  ReadProcpar('gror',filename_procpar);  % Read out prephasing in G/cm
% procpar.tpe =  ReadProcpar('tpe',filename_procpar); % Phase encode gradient time in seconds
% procpar.np = ReadProcpar('np',filename_procpar); % Number of time points along readout direction (real and imaginary)
% procpar.nv = ReadProcpar('nv',filename_procpar); % Number of phase encode lines
% procpar.orient = ReadProcpar('orient',filename_procpar);
% procpar.petable_name = ReadProcpar('petable',filename_procpar);
% if strcmp(procpar.petable_name,'n') == 0
%     procpar.petable = linspace(procpar.nv/2,-(procpar.nv/2)+1,procpar.nv);
% end

%%

load procpar.mat

% Bfield and phantom------------------------------------------------------
m = phantom(32,32);
n=procpar.nv;

%Original field map
load 'C:\Users\Vic\Documents\MATLAB\LabWork\TMS_AC\wire_circ_32.mat';

%Define the static magnetic field B0
B0x = zeros(n,n);
B0y = zeros(n,n);
B0z = ones(n,n) .* 1.5 .* 1e4;


% biot_3D does calculations in Tesla. All calculation in this code are done
% in gauss. 1 T = 10^4 gauss
Bfx = bx(:,:,12) .* 1e4;
Bfy = by(:,:,12) .* 1e4;
Bfz = bz(:,:,12) .* 1e4;

% Bfield = abs(B0 + B2) - abs(B0), where B2= DC field map
babs = sqrt( (Bfx + B0x) .^ 2 + (Bfy + B0y) .^ 2 + (Bfz + B0z) .^ 2 );
B0_mag = sqrt( B0x .^ 2 + B0y .^ 2 + B0z .^ 2 );

amp = 10; % Manipulate amplitude according to requirements, If amp=1 current = 1mA
current = amp*1e-3

Bfield_amp = (babs - B0_mag ) .* amp;


% fovxy = 20; % field of view
% xmax = fovxy/2;
% ymax = fovxy/2;
% xstep = fovxy/n;
% ystep = fovxy/n;
% [x y] = meshgrid(-xmax:xstep: xmax-xstep  , -ymax:ystep: ymax-ystep );
% Bfield_amp = 0.001 * exp(-((x).^2 + (y+3).^2)/5);




clear B0x BOy B0z Bfx Bfy Bfz babs B0_mag bx by bz
%--------------------------------------------------------------------------
%%
[kro,kpe,xx_mat,yy_mat] = kspace_gems(procpar);
procpar.xres=procpar.lpe/procpar.nv;
procpar.yres=procpar.lro/(procpar.np/2);

procpar.freq = 1e3;
procpar.initial_phase = 0;
[phase_vec,swave] = phase_calc_GE(procpar);

procpar.kro = kro;
procpar.kpe = kpe;
procpar.phase_vec = phase_vec';
procpar.xx_mat = xx_mat;
procpar.yy_mat = yy_mat;
procpar.m = m;

%%
% Calculation for signal without Bfield----------------------------------
index = find(abs(m)>0);
procpar.index = index;
procpar.sin_phase = pi/2;
Bfield_r = zeros(size(index));

signal_field_off = signalDC_lsq(Bfield_r,procpar);


tmp = reshape(signal_field_off,n,n);
image_field_off = ift2(tmp); %% Use this as procpar.m

%%
% Calculation for signal with Bfield

% procpar.weight = ones(procpar.nv,procpar.nv);
Bfield_r =Bfield_amp(index) ;
procpar.sin_phase = 0;
signal_field_on1 = signalDC_lsq(Bfield_r,procpar);

tmp = reshape(signal_field_on1,n,n);
% tmp = fliplr(flipud(tmp));
image_field_on1 = ift2(tmp);


%%
% Calculation for signal with Bfield

Bfield_r =Bfield_amp(index) ;
procpar.sin_phase = pi/8;
signal_field_on2 = signalDC_lsq(Bfield_r,procpar);

tmp = reshape(signal_field_on2,n,n);
% tmp = fliplr(flipud(tmp));
image_field_on2 = ift2(tmp);


%%

phase1 = angle(conj(image_field_on2) .* image_field_on1);
Bf_est = zeros(n,n);
Bf_est(index) = phase1(index);

Bf_est = Bf_est ./ (2 * pi * cos(procpar.sin_phase) - cos(0));
figure;
subplot(121),imagesc(Bf_est);
subplot(122),imagesc(Bfield_amp);





