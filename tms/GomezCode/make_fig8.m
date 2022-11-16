function [wire Jvec] = make_fig8(current, r, zpos, theta)
%function make_fig8(current, r, zposition, bend_theta)


    FOV = 0.5; % m
    

    wire =[];
    
    % Make the fig8 coils:
    loopSegs = linspace(0,2*pi,20);
   
    wx = r*sin(loopSegs);
    wy = -r*cos(loopSegs);
    wz = zpos + wy*sin(theta);
   
    loop =[wx' -r-wy' wz';
        wx' r+wy' wz'];

    wire = [wire; loop];
    l=length(wire(:,1));
    Jvec(l,3)=0;
    for i=1:l
        Jvec(i,:)=wire(1+i*(i~=l),:)-wire(i,:);
    end
return