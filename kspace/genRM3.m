function R = genRM3(axises,angles)
%function genRM3(axises,angles)
%|
%|  Generates a 3x3 rotation matrix given axis/angles of euler rotations
%|
%|  Inputs:
%|      axises: array of strings with given axises in multiplication order
%|          (example: ["y" "x" "z"] for Ry*Rx*Rz)
%|      angles: array of angles in multiplication order
%|          (example: [pi/2 0 pi])
%|

    % Initialize R as an identity matrix:
    R = eye(3);

    % Specific axis rotation matrices:
    Rx = @(xi) [1 0 0; 0 cos(xi) -sin(xi); 0 sin(xi) cos(xi)];
    Ry = @(psi) [cos(psi) 0 sin(psi); 0 1 0; -sin(psi) 0 cos(psi)];
    Rz = @(phi) [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

    % Loop through inputs and create a multiplicative sum of matrices:
    for i = length(axises):-1:1
        axis = axises(i);
        angle = angles(i);
        switch axis
            case 'x'
                R = Rx(angle)*R;
            case 'y'
                R = Ry(angle)*R;
            case 'z'
                R = Rz(angle)*R;
        end
    end
end

