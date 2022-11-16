function [Gmap, psi] = SENSE_simulation
%
% (c) 2009 - Luis Hernandez-Garcia @ UM
% Last Edit Dec. 15, 2009
% 
% do a simulated SENSE reconstruction on cartesian data


% 
% allS = read_BFields('metaSims128Mhz/coilnometa',8);
% %allS = allS(:, 1:2:end); % downsample to 64x128
% psi = coupling16rungs;
% 
% Gmap = makeGmap(allS, psi);

N=32;
FF = FourierMat(N);  % 2D FT matrix

% close all
 R = 2;
Fdown = zeros((N^2)/R, N^2);
for n=1:N/R
    start1= (n-1)*N +1; fin1 = (n-1)*N +N;
    start2= (n-1)*N*R +1;  fin2 = (n-1)*N*R +N;
    Fdown(start1 :fin1, :) = FF( start2 : fin2, :);
end
imagesc(real(Fdown));

x=phantom(N);
imagesc(x);
x= x(:);
xf = FF*x;
imagesc(abs(reshape(xf,N,N)))

x2 = Fdown* x; 
x2 = reshape(x2, N, N/2);
imagesc(abs(x2));

% zero fill in k-space:
x3 = [zeros(N,N/4)  x2  zeros(N,N/4)];
imagesc(abs(x3));
x3 = x3';

%FT again:
x4 = F*x3(:);
imagesc(reshape(abs(x4),32,32));


keyboard



%%
function FF=FourierMat(N)
% function F=FourierMat(N)
%
% builds a complex Fourier coefficient matrix of size N
% you can do and FFT of x by doing:  F*x
%

w = linspace(-pi, pi-2*pi/N, N);
F = zeros(N);
t = [0:N-1];
for n=1:N
    F(n,:) = exp(-i*w(n)*t)/(N*N); 
end

% turns out matlab already had that built in!
% F = dftmtx(N);

FF=zeros(N*N);  % now the 2D FT matrix
for n=1:N
      FF( (n-1)*N+1: n*N, (n-1)*N+1: n*N) = F * exp(-i*w(n));
end

%FF=repmat(F,N,N);
FF=zeros(N*N);
for n=1:N
    for m=1:N
        FF( (m-1)*N+1: m*N, (n-1)*N+1: n*N) = F * exp(-i*w(m)*w(n)*pi);
    end
end

imagesc(real(FF))

x=phantom(N);
imagesc(x);
x= x(:);
xf = FF*x;
imagesc(abs(reshape(xf,N,N)))
xff = FF*xf;
imagesc(abs(reshape(xff,N,N)))

return


%%
function [allS ] = read_BFields(basestr, N)
% this is the function that reads in the .mat files with the simulated B
% fields

allS = [];
for n=1:N
    str = [basestr num2str(n) '.mat']
    
    load(str);
    Bx = abs(Bx); 
    By = abs(By); 
    B = complex(Bx,By);
    figure(1); imagesc(abs(By)); drawnow, pause(0.1);
    B = B(1:4:end, 1:4:end);
    allS = [allS; B(:)'];
end


%%
function Gmap = makeGmap(allS, psi)
% function Gmap = makeGmap(allS, psi)
%
% this computes the G factor map from sensitivitiy maps and noise
% correlation matrix from Pruessmann:  Magnetic Resonance in Medicine 42:952?962 (1999)
% 
%      Gmap(:) = diag(sqrt( inv(SH * psi_inv * S) .* (SH * psi_inv * S)));

Ncoils = size(allS,1);
Npix = size(allS,2);

psi_inv = pinv(psi);
psi_inv=eye(8);

Gmap = zeros(1,Npix);

% first try
% for p=1:Npix
%     S = allS(:,p);
%     SH = transpose(conj(S));
% 
%     Gmap(p) = ( inv(SH * psi_inv * S) .* (SH * psi_inv * S)).^0.5 ;
% end

% second try:
allSH = transpose(conj(allS));
Gmap = (diag( inv(allSH * psi_inv * allS) .* (allSH * psi_inv * allS))).^0.5 ;


Gmap = reshape(Gmap, sqrt(Npix), sqrt(Npix));
figure(3); imagesc(abs(Gmap));


return




%%
function couplingMat = coupling16rungs
% function couplingMat = coupling16rungs
% 
% this function computes the inductance of 8 infinite square loops around a
% cylinder.  The rungs (16 of them) are equidistant and the current
% alternates between them
% This is done in 2D - a cross section through an infinite cyllinder
%
% we compute it by calculating the B field produced by each coil (source) at the
% surface of the other coils (targets), then we use Faraday's induction law by
% integrating that field on the surface and computing the induced EMF on
% each target coil
%
% the result is a matrix of EMFs at each coil induced by every other coil.
% The result is normalized by the diagonal.
%

dIdt=1;
R=0.15;
coilx = linspace(0,2*pi - 2*pi/16, 16);
coilx = R*sin(coilx);
coilx = [coilx coilx(1)];

coily = linspace(0,2*pi - 2*pi/16, 16);
coily = R*cos(coily);
coily = [coily coily(1)];

%plot(coilx, coily);


couplingMat = zeros(8,8);
for n= 1:2:16
    
    if n==1, doShow=0, else doShow=0, end;
    
    source = [
        coilx(n) coily(n) 0;
        coilx(n+1) coily(n+1) 0
        ];
        %dIdt = dIdt * (-1)^n;

    for m=1:2:16
        target = [
            coilx(m) coily(m) 0;
            coilx(m+1) coily(m+1) 0
            ];
        couplingMat(round(n/2),round(m/2)) = induction2D(dIdt,source,target, doShow);
    end
    
end

%for n=1:8, couplingMat(n,n)=0; end
couplingMat = couplingMat/couplingMat(1,1);

figure(2)
imagesc(log10(couplingMat)); drawnow;
return



%%
function volts = induction2D(dI_dt, source, target, showVectors)
% the units have to be in meteres and amps/sec.

% number of places where to sample the magnetic field
Nelements = 10;
x = linspace(target(1,1) , target(end,1), Nelements);
y = linspace(target(1,2) , target(end,2), Nelements);
z = zeros(size(x));
locations = [x; y ; z]';

MU = 1.26e-6;  % Permeability H/m
% this is the term that goes into Faraday's Law
K = -dI_dt * MU / (4*pi);
B = [0 0 0];

for r=1:Nelements

    % this is the vector from the source rung to the point in the coil
    r1 = source(1,:) - locations(r,:) ;
    % the current is along the z-axis
    dl1 = [0,0,1];
    D1 = sum((locations(r,:) - source(1,:)).^2)^0.5;
    B1 = cross(r1 , dl1) / (D1+eps)^2;


    % do the same thing for the second rung, whose current runs in the
    % opposite direction
    r2 = source(2,:) - locations(r,:) ;
    dl2 = [0,0,-1];
    D2 = sum((locations(r,:) - source(2,:)).^2)^0.5;
    B2 = cross(r2 , dl2) / (D2+eps)^2;

    if showVectors
        % option: visualize the B-field vectors induced by each rung:
        axis (0.3*[-1 1 -1 1])
        quiver(x(r), y(r), B1(1)/100, B1(2)/100,'r'); hold on; quiver(x(r), y(r), B2(1)/100, B2(2)/100,'g');  drawnow
    end
    
    % integrate the fields
    B = B + B1 + B2;

end

% this is the integral of the field on that surface (or line in this case)
B = B/Nelements;

% Compute the normal component of this vector relative to the
% target coil's plane?
targetplane = target(2,:) - target(1,:);
Bmag = norm(B);
Bproj = B * targetplane';
Bnormal = Bmag*sin(acos(Bproj/Bmag));

% use Faraday's law to scale it appopriately
volts = K * Bnormal;

return
