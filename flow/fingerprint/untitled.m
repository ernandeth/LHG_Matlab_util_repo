% Radial Sampling script :  how to choose angles for radial sampling
% rotations in 3D

N = 60;
p = [1,0, 0];
doFib=0

angz = deg2rad(137.5);
angy = 0;
angx = 3*pi/N;  % option 1

GRatio = (1+sqrt(5))/2;

%angx = angz;  % option 2

randos = rand(100,1);
close

for n=0:N-1

    if doFib
        angx = 2*pi*n/GRatio;
        angz =  acos(1-2*(n+0.5)/N);
    end

    rotz=[cos(angz*n) -sin(n*angz) 0; 
        sin(n*angz) cos(n*angz) 0 ; 
        0 0 1];
    
    rotx=[1 0 0; 
        0 cos(angx*n) -sin(n*angx); 
        0 sin(n*angx) cos(n*angx)];
    
    roty=[ cos(angy*n) 0 sin(n*angy) ; 
        0 1 0;  
        -sin(n*angy) 0 cos(n*angy) ];
    
    pp = roty*rotx*rotz*p';

    line([pp(1) -pp(1)], [pp(2) -pp(2)], [pp(3) -pp(3)]);
    
    norm(pp);
    
    hold on
    drawnow

end
xlabel('x')
ylabel('y')
zlabel('z')
