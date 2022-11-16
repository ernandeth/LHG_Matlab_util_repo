%   Function to rotate k space trajectoies using a simple rotation matrix
%   at angles spaced equally about the unit circle: [0 pi/Nint,
%   2*pi/Nint, ..., (Nint-1)*2*pi/Nint]
%   [trajVec, trajMat] = rotTraj(Trajectory, nInt);
%        Inputs:   Trajectory: is a nx2 vector: [xcoords,ycoords];
%                  Nint:  The number of interleaves used during imaging
%        Outputs:  trajVec: A (length(phaseAngles)*n)x2 vector where each 
%                           rotated trajectory is concatenated onto trajVec
%                   trajMat:  A nXlength(phaseAngles)X2 matrix containing
%                             the rotated trajectories in the columns. 
%                             EX: to access xcoords rotated by 2*Pi*n/Nint radians 
%                             enter xTrajn = trajMat((:,n,1);
function    [trajVec, trajMat] = rotTraj(Trajectory, Nint)

% Set up constants             
[m n] = size(Trajectory);

% Ensure that the vector is nx2 for rotation
if n>m;
    Trajectory = Trajectory';
    [m n] = size(Trajectory);
end


% Preallocate Memory
trajMat = zeros(m,Nint,2);  % Preallocate memory
trajMat(:,1,1) = Trajectory(:,1);
trajMat(:,1,2) = Trajectory(:,2);

% Rotate the trajectories
phase = (0:Nint-1)*2*pi/Nint;

for q = 2:size(phase,2)

%  Set up Rotation Matrix.
  tmpTraj = rotclkwise(Trajectory, -phase(q));  % Perform Rotation (nx2x2x2)
%   tmpTraj = Trajectory*rotclkwise;  
  trajMat(:,q,1) = tmpTraj(:,1);  %  X coordinate Storage
  trajMat(:,q,2) = tmpTraj(:,2);  %  Y coordinate Storage
end

% Create trajVec
trajVec = [reshape(trajMat(:,:,1),m*Nint,1), reshape(trajMat(:,:,2),m*Nint,1)];
end
                   
