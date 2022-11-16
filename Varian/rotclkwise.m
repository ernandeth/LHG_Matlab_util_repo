% %   Function to rotate a vector counter-clockwise with respect to the x,y coordinate
% frame.
%      function returnVector = rotclkwise(inputVector, angle)
%      Inputs:   inputVector:  a nx2 vector in the form [x coords, y
%      coords].
%                angle:  angle of rotation in radians
%      Outputs:  returnVector = inputVector*rotationMatrix(angle);
function   returnVector = rotclkwise(inputVector, angle)
              
%   make sure vector is nx2
    [m,n] = size(inputVector);
    if n>m
        inputVector = inputVector';
        display('inputVector should be nx2')
        display('    And in the form [x coords, ycoords]')
    end
  rotationMatrix = [cos(angle), -sin(angle); sin(angle),cos(angle)]'; % Clockwise rotation
  returnVector = inputVector*rotationMatrix';