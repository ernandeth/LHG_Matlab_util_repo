function [Es Esx Esy Esz] = make_square_shield(current, R, zpos) 
%function [Es Esx Esy Esz] = make_square_shield(current, R, zpos) 
    
global FOV Nvox

doSwapxy=0;
wire=[];



wy = [0;
    
    linspace( R/20, R,  5)';
    R*[1 1 1 1 1]';
    linspace(R, R/20, 5)';
    
    0
    
    linspace( -R/20, -R,  5)';
    -R*[1 1 1 1 1]';
    linspace(-R, -R/20, 5)';
    
    0
    ];

wx = [0;
    
    -R/2*[1 1 1 1 1]';
    linspace( -R/2, R/2,  5)';    
    R/2*[1 1 1 1 1]';
    
    0;
    
    -R/2*[1 1 1 1 1]';
    linspace( -R/2, R/2,  5)';    
    R/2*[1 1 1 1 1]';
    
    0
    ];

wz =  zpos* ones(size(wx));
wz([1,17,33])=zpos+10;
% 
% for count=1:length(wx)
%     plot3(wx(count), wy(count), wz(count),'*'); hold on
%     axis (1.1*[-R R -R R zpos-1 zpos+1])
%     pause(0.1)
% end

wire =[wx wy wz];

%%%
if doSwapxy
    tmp=wire;
    tmp(:,1) = wire(:,2);
    tmp(:,2) = wire(:,1);
    wire=tmp;
end
%%
sfigure(3);
plot3(wire(:,1), wire(:,2), wire(:,3),'b');
axis([-0.25 0.25 -0.25 0.25 -0 0.25]);
hold on
drawnow


% Calculate E field from coil alone:
[Es Esx Esy Esz] =  Efield(current, wire, [Nvox Nvox Nvox], FOV);



return



