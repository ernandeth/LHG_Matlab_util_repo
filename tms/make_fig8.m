function [E1 E1x E1y E1z] = make_fig8(current, r, xgap, zpos, theta)
%function make_fig8(current, r, zposition, bend_theta)

global FOV Nvox


    wire =[];
    
    % Make the fig8 coils:
    loopSegs = linspace(0,2*pi,10);
   
    wx = r*sin(loopSegs) ;
    wy = -r*cos(loopSegs);
    wz = zpos - wy*sin(theta) - r*sin(theta) ;
   
    loop =[wx'  -r-wy'-xgap/2 wz';
	0 0 10;
        wx' r+wy'+xgap/2 wz'];

    wire = [wire; loop];

	%sfigure(3);
    %plot3(wire(:,1), wire(:,2), wire(:,3),'k')
%     hold on , plot3(0,0,0,'ro')
%     axis([-0.25 0.25 -0.25 0.25 -0 0.25]);
%     hold on
%     drawnow
    % Allocate space and define some conductivity space
    
    % Calculate E field from coil alone:
   [E1 E1x E1y E1z] =  Efield(current, wire, floor([Nvox Nvox Nvox]),FOV);
return
