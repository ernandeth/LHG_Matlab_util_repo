%   Function to reconstruct spiral data

%function  reconstructSpirals(studyPath,spiralScanName, KxTrajectory, KxReference,KyTrajectory, KyReference);

function  picture = reconstructSpirals(varargin)  
%reconstructSpirals(studyPath,spiralScanName, KxTrajectory,KyTrajectory);



% Parse inputs
td = varargin{1};               % Top directory
spiralData = varargin{2};
reps = varargin{3};
if numel(varargin)>=4
KxDir = varargin{4};
KyDir = varargin{5};
end

% Set Constants
gamma = 4257;  %radians/s/G
dt = 4e-6;     % Gradient sample time in seconds

% Get spiral data.
[dta, procpar] = getSpiralFidEpi([td '/' spiralData], reps);
sw = procpar.sw;
ts = (1:size(dta,1))/sw;

% get Kspace locations and upsample to correct number of points if
% necessary
if numel(varargin) >=4;
KyMeas = getTraj([td '/' KyDir]);
KxMeas = getTraj([td '/' KxDir]);
end
[Gyi, Kyi, tg] = getIdealGrads([td '/' spiralData '/Waveforms'], 'Gy.GRD', dt, sw, 1:size(dta,1));
[Gxi, Kxi, tg] = getIdealGrads([td '/' spiralData '/Waveforms'], 'Gx.GRD', dt, sw, 1:size(dta,1));

Kxi  = Kxi - Kxi(end);
Kyi  = Kyi - Kyi(end);

% Rotate trajectories for multishot imaging
[Ki, ~] = rotTraj([Kxi,Kyi],procpar.Nint);
[Gi, ~] = rotTraj([Gxi,Gyi],procpar.Nint);

if numel(varargin)>=4 
[KMeas,~] = rotTraj([KxMeas, KyMeas],procpar.Nint);
end


% Regrid Data
Regrid = zeros(2*procpar.Mat, 2*procpar.Mat, procpar.ns, reps);  % Preallocate Memory
for l = 1:reps
for q = 1:procpar.ns
    if numel(varargin) < 4
    Regrid(:,:,q,l) = matlabRegrid(Ki,dta(:,l,q),procpar.lro, procpar.Mat,procpar.Mat*2, procpar);
    end
    if numel(varargin) >= 4
    Regrid(:,:,q,l) = matlabRegrid(KMeas,dta(:,l,q),procpar.lro, procpar.Mat,procpar.Mat*2, procpar);
    end
    % Compute the ifft
    picture(:,:,q,l) =ifftshift(ifft2(Regrid(:,:,q,l))) ;
end
l
end

end




