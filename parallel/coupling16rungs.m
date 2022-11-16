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
    
    if n==1, doShow=0; else doShow=0; end;
    
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
