function out = ellipse_mask(in, rx, ry, rz)
% function out = ellipse_mask(in, rx, ry, rz)
% 
% mask the data of a 3D volume such that everything in the
% keep everything inside an ellipse determined by rx, ry, rz
% and centered on the center of gravity of the image
%
% returns the mask (- not the data)
%

out = in;
out(:) = 0;

[xdim ydim zdim] = size(in);
Npix = xdim*ydim*zdim;

% find center of mass of the data:
xcenter = 0;
ycenter = 0;
zcenter = 0;
total_mass = sum(in(:));
for x=1:xdim
    for y=1:ydim
        for z=1:zdim
            weight = in(x,y,z)/total_mass;
            xcenter = xcenter + x*weight;
            ycenter = ycenter + y*weight;
            zcenter = zcenter + z*weight;
        end
    end
end

% center of mass:
com = [xcenter ycenter zcenter];

for x=1:xdim
    for y=1:ydim
        for z=1:zdim
            if ( ((x-xcenter)/rx)^2 + ((y-ycenter)/ry)^2 + ((z-zcenter)/rz)^2 <1 )
                out(x,y,z)=1;
            end
        end
    end
end

return
