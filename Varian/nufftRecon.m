%Function to reconstruct spiral data using Jeff Fessler's IRT programs.
%Uses the fatrix2 formulation for iNUFFT regridding and weighted least
%squares with roughness penalization.
%Installation of IRT is prerequisit before using this program.

%nufftRecon(studyPath,spiralScanName, KxTrajectory,KyTrajectory);
%   Inputs: studyPath: String. Top directory of study folder containg the .fid folders of the scans in question.
%                       Ex:  '/home/<user>/vnmrsys/data/<study name>';
%           sprialScanName: String. Name of the .fid folder of the scan you wish to reconstruct. 
%                           Ex: 'spiral_01.fid';
%           KxTrajectory: (optional) String. Name of the .fid folder containing Kx encoding FID data.
%                         Ex: 'KxMap_02.fid';
%           KyTrajectory: (optional) String. Name of the .fid folder containing Ky encoding FID data.
%                         Ex: 'KyMap_02.fid';
%   Outputs: picture: Inverse Fast Fourier Transform of Regridded Kspace data.

function  picture = nufftRecon(varargin)  
%  nufftRecon(studyPath,spiralScanName, KxTrajectory,KyTrajectory);



% Parse inputs
td = varargin{1};               % Top directory
spiralData = varargin{2};       % file holding the acquired data
if numel(varargin)>=3
KxDir = varargin{3};            % Directory holding the mesured and refence trajectories.
KyDir = varargin{4};
end

% Set Constants
gamma = 4257;  %radians/s/G
dt = 4e-6;     % Gradient sample time in seconds

% Get spiral data.
[dta, procpar] = getFid([td '/' spiralData]);
dta2 = reshape(dta,numel(dta),1);
sw = procpar.sw;
ts = (1:size(dta,1))/sw;

% get Kspace locations and upsample to correct number of points if
% necessary
[Gyi, Kyi, tg] = getIdealGrads([td '/' spiralData '/Waveforms'], 'Gy.GRD', dt, sw, 1:size(dta,1));
[Gxi, Kxi, tg] = getIdealGrads([td '/' spiralData '/Waveforms'], 'Gx.GRD', dt, sw, 1:size(dta,1));

if numel(varargin) >=3;
KyMeas = getTraj([td '/' KyDir]);
KxMeas = getTraj([td '/' KxDir]);
end

% Rotate trajectories for multishot imaging
[Ki, ~] = rotTraj([Kxi,Kyi],procpar.Nint);
[Gi, ~] = rotTraj([Gxi,Gyi],procpar.Nint);

if numel(varargin)>=3  
[KMeas,~] = rotTraj([KxMeas, KyMeas],procpar.Nint);
end

% Setup Nufft reconstruction
kspace = KMeas;
ig = image_geom('nx', procpar.Mat, 'ny', procpar.Mat, 'dx', procpar.lro/procpar.Mat, 'offsets', 'dsp');
N = ig.dim;
J = [6 6];
[a b c] = reconstructGREImage([td '/gems_10.fid']);
U = dtft2(abs(a), 2*pi*KMeas*c.lro/c.nv, [16 16]);
ts = 1:length(U);
plot(ts,abs(U)/max(U(:)), ts, abs(KMeas(:,1))/max(KMeas(:,1)), ts, abs(KMeas(:,2))/max(KMeas(:,2)));

nufft_args = {N, J, 2*N, N/2, 'table', 2^10, 'minmax:kb'};
Gm = Gmri(kspace, ig.mask, ...
		'fov', ig.fov, 'basis', {'rect'}, 'nufft', nufft_args);


%  Reconstruction
% 	Rn = Robject(ig.mask, 'type_denom', 'matlab', ...
% 		'potential', 'hyper3', 'beta', 2, 'delta', 0.3);
%     xh = Gm'*(dta2);   % Conjugate phase reconstruction.
	xh = qpwls_pcg1(Gm'*dta2, Gm, 1, dta2(:), 1.5, 'niter', 120);

    recxh = reshape(xh,[ig.dim]);
    
% Account for lateral translations of the desired Field of View
Regrid = fftshift(fft2(recxh));
Mat = procpar.Mat;
FOV = procpar.lro;
gridSize = Mat;
x0 = procpar.pro;
y0 = procpar.ppe;
KXi = linspace(-Mat/FOV/2,Mat/FOV/2,gridSize);
[KXi,KYi] = meshgrid(KXi,KXi);
fx = 1i*2*pi*x0*KXi;
fy = 1i*2*pi*y0*KYi;
Regrid1 = Regrid.*exp(fx+fy).*exp(fy);

picture = (ifft2(Regrid1));
    
end


