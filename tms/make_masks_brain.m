% make a headsphere
    headmask = zeros(size(E1));
    
    dims = size(headmask)
    center = size(E1)/2
    head_radius = min(center)-5;
    
    for x=1:dims(1), for y=1:dims(2), for z=1:dims(3)
                if sqrt((x-center(1))^2 +(y-center(2))^2 + (z-center(3))^2) <= head_radius,
                    headmask(x,y,z) = 1;
                end,
            end,end,end


    % make a target
    targetmask = zeros(size(E1));
    center(3) = 20;
    target_radius = 2;
    for x=1:dims(1), for y=1:dims(2), for z=1:dims(3)
                if sqrt((x-center(1))^2 +(y-center(2))^2 + (z-center(3))^2) <= target_radius,
                    targetmask(x,y,z) = 1;
                end,
            end,end,end


    % make a cylinder
    cylmask = zeros(size(E1));
    %head_radius = 16;
    
    for x=1:dims(1), for y=1:dims(2), for z=1:dims(3)            
                if sqrt((x-center(1))^2 +(y-center(2))^2 ) <= 1,
                    cylmask(x,y,z) = 1;
                end,
            end,end,end
    cylmask = cylmask .* headmask;
