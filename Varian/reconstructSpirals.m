%   Function to reconstruct spiral data

%function  reconstructSpirals(studyPath,spiralScanName, KxTrajectory, KxReference,KyTrajectory, KyReference);

function  picture = reconSpOutInEngine(varargin)  
%reconstructSpirals(studyPath,spiralScanName, KRoTrajectory,KPeTrajectory);



% Parse inputs
td = varargin{1};               % Top directory
spiralData = varargin{2};
if numel(varargin)>=3
KRoDir = varargin{3};
KPeDir = varargin{4};
end

% Set Constants
gamma = 4257;  %radians/s/G
dt = 4e-6;     % Gradient sample time in seconds

% Get spiral data.
[dta, procpar] = getSpiralFid([td '/' spiralData]);
% dta = dta./max(abs(dta));
sw = procpar.sw;
ts = (1:size(dta,1))/sw;

% get Kspace locations and upsample to correct number of points if
% necessary
if numel(varargin) >=3;
KRoMeas = getTraj([td '/' KRoDir]);
KPeMeas = getTraj([td '/' KPeDir]);


end
[Gpei, Kpei, Groi, Kroi, ts, GPe, KPe, GRo, KRo, tg] = getIdealGrads(procpar);

% Rotate trajectories for multishot imaging
[Ki, ~] = rotTraj([Kpei,Kroi],procpar.Nint);
[Gi, ~] = rotTraj([Gpei,Groi],procpar.Nint);

if numel(varargin)>=3  
[KMeas,~] = rotTraj([KPeMeas, KRoMeas],procpar.Nint);
end

dta2 = dta;

% Regrid Data
Regrid = zeros(2*procpar.Mat, 2*procpar.Mat, procpar.ns);  % Preallocate Memory
picture = zeros(size(Regrid));
for q = 1:procpar.ns
    if numel(varargin) < 3
    Regrid(:,:,q) = matlabRegrid(Ki,(dta2(:,:,q)),procpar.lro, procpar.Mat,procpar.Mat*2, procpar);
    end
    if numel(varargin) >= 3
    Regrid(:,:,q) = matlabRegrid(KMeas,(dta2(:,:,q)),procpar.lro, procpar.Mat,procpar.Mat*2, procpar);
    end
    picture(:,:,q) = ifftshift(ifft2((Regrid(:,:,q))));  % Compute the ifft

end



end




