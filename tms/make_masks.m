global Nvox, FOV;

% make a headsphere
headmask = zeros(Nvox, Nvox, Nvox);
headmasksurf = zeros(Nvox, Nvox, Nvox);

zshift=0.02*Nvox/FOV;
center = ceil([Nvox/2, Nvox/2, (Nvox/2)]);
%center = ceil([Nvox/2, Nvox/2, Nvox/2-3]);  % this is 4.53 cm down from "isocenter"


head_radius = 0.075;  %meters
head_radius = floor(head_radius*Nvox/FOV); % pixels

for x=1:Nvox, for y=1:Nvox, for z=1:Nvox
            if abs(head_radius - sqrt((x-center(1))^2 +(y-center(2))^2 + (z-center(3))^2) )<= 0.5,
                headmasksurf(x,y,z) = 1;
            end,
            if sqrt((x-center(1))^2 +(y-center(2))^2 + (z-center(3))^2) <= head_radius,
                headmask(x,y,z) = 1;

            end,
        end,end,end





% make a cylinder
cylmask = zeros(size(headmask));
% head_radius = Nvox/4;
% center = round([Nvox/2, Nvox/2, Nvox/2-4]);
for x=1:Nvox, for y=1:Nvox, for z=center(3)-head_radius :center(3)+head_radius
            z=round(z);
            
            % this is in pixels
            if sqrt((x-center(1))^2 +(y-center(2))^2 ) <= 0.5,
                cylmask(x,y,z) = 1;
            end,
        end,end,end


% make a target
targetmask = zeros(size(headmask));
center = [Nvox/2, Nvox/2, Nvox/2];
target_radius = 1;
for x=1:Nvox, for y=1:Nvox, for z=1:Nvox
            if sqrt((x-center(1))^2 +(y-center(2))^2 + (z-center(3))^2) <= target_radius,
                targetmask(x,y,z) = 1;
            end,
        end,end,end

if 0
    ov([], (1+headmask+cylmask)*10 ,...
        ceil(Nvox/2), ceil(Nvox/2), ceil(Nvox/2),0);
end
