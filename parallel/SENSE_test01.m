% function [Gmap, psi] = SENSE_simulation
%
% (c) 2009 - Luis Hernandez-Garcia @ UM
% Last Edit Dec. 15, 2009
%
% do a simulated SENSE reconstruction on cartesian data
% this works for simplistic downsampling


N=128;
Npix = N^2;
Ncoils = 8;
x = zeros(N);
x(N/4+1:end-N/4, N/4+1:end-N/4) = phantom(N/2);
%x(1:10,N-10:N)=10;
x = x(:);

if 1
figure
subplot(221)
allS = read_BFields('metaSims128Mhz/coilmeta',N);
% allS = fake_Bfields(N);

psi = coupling16rungs;
psi = eye(8);
end

R = 8;
Gmap = zeros(1,Npix);
x_sense = zeros(1,Npix);
alias_gap = ceil(N/R)+1;
ytmp = zeros(Ncoils, Npix);
warning off

for col=1:N
    for row=floor((N-N/R)/2) : floor((N +N/R)/2)
        % figure out which pixels contribute to alias and extract those
        % sensitivities.  (aliasing occurs in ydirection)

        p = sub2ind([N, N], row, col);
        row_alias = [];
        
        for r=-R:R
            row_alias =  [row_alias ; row + r*alias_gap];
        end

        % make sure I'm in the FOV
        row_alias = row_alias(row_alias>0);
        row_alias = row_alias(row_alias<N);
        col_alias = col*ones(size(row_alias));

        % now get their indices in the array - we note that rows = y-coordinate
        p_alias = sub2ind([N, N], row_alias, col_alias);
        
        % encoding goes like this:  y = S * x
        % need to decode by x = inv(S)*y

        % sensitivities at locations that alias on top of each other
        Stmp = allS(:,p_alias);  
        % the aliasing pattern going forward
        ytmp(:,p) = Stmp * x(p_alias);
        
        % try to recover x by inversion of Stmp
        % xtmp = pinv(Stmp)*ytmp(:,p);

        U = makeUnfoldMat(Stmp, psi);
        xtmp = U*ytmp(:,p);
        
        %xtmp(isnan(xtmp)) = 0;
        %xtmp(isinf(xtmp)) = 0;
        
        x_sense(p_alias) =  xtmp' + x_sense(p_alias) ;

        Gmap(p_alias) = makeGmap(allS(:,p_alias), 1);
    end
end

figure
subplot(222)
for n=1:Ncoils
    imagesc(reshape(abs(ytmp(n,:)),N,N));    
    title(sprintf('aliased image from coil %d', n));
    drawnow, pause(0.1)
end

subplot(223)
imagesc(reshape(abs(x_sense), N,N));
title('Reconstructed image')

subplot(224)
imagesc(reshape((Gmap), N,N)); 
colorbar('EastOutside')
title('G factor Map')

