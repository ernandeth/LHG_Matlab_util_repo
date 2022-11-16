%   Function to get Measure Kx, Ky trajectories for a particular object.
%  [KTraj] = getTraj(fileLocation)
%    Inputs: filteLocation: String. Location of the .fid file produced by Vnmrj
%                           which contains the FID of the desired scan.
%                           Ex:/Home/<User>/vnmrsys/data/<study>/<scan>.fid.
%    Trajectories are obtained by subtracting a the phase of the signal
%    from two slices seperated along one dimension.  See
%    (Duyn JG, et al. May 1998. Simple Correction Method for K-Space Trajectory
%     Deviations in MRI. J Magn Reson. 132(1):150:3.)

function [KTraj] = getTraj(fileLocation)
% [KTraj] = getTraj(fileLocation, referenceLocation, dgap, sw)

% Set up De noising filter for trajectory data.  (We assume K-space trajectory is
% smooth and continuous)
% gfilter = gausswin(15);
% gfilter  = gfilter./(sqrt(gfilter'*gfilter));  % Normalize

% get the trajectory data
[u, procpar] = getFid(fileLocation);
[a b] = size(u);
idx1 = 1:b/2;
idx2 = b/2+1:b;

% Compute the trajectory for each slice
% if strcmp(procpar.spiral_in,'y')
%     u = u(end:-1:1,:);
% end


% if spiral in, unwrap from the back in
% if strcmp(procpar.spiral_in,'y')
%     u = u(end:-1:1,:);
% end

% now unwrap the phase
kref = unwrap(angle(u(:,idx2)),pi,1);
kRaw = unwrap(angle(u(:,idx1)),pi,1);

k = kRaw-kref;

% if spiral in, restore the direction of the waveform
% if strcmp(procpar.spiral_in,'y')
% k = k(end:-1:1,:);
% end

% calculate trajectsory
scalingMatrix = diag(1./procpar.pss)/2/pi;
KTraj = k*scalingMatrix;

% if strcmp(procpar.spiral_in,'y')
% KTraj = KTraj-ones(size(KTraj))*diag(KTraj(end,:));
% end

% Introduce some way to throw out outliers where the trajectory is just too
% off to use.

idealMax = procpar.Mat/procpar.lro/2;
if strcmp(procpar.spiral_in,'y')
    ratioMax = max(KTraj(200:end,:))/idealMax;
    ratioMin = min(KTraj(200:end,:))/(-idealMax);
else
    ratioMax = max(KTraj(1:end-200,:))/idealMax;
    ratioMin = min(KTraj(1:end-200,:))/(-idealMax);
end

idxs = find(ratioMax >= 0.7 & ratioMax <= 1.3 & ratioMin >= 0.7 & ratioMin <= 1.3);


KTraj = mean(KTraj(:,idxs),2);



end
