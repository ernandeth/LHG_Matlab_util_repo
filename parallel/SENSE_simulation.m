% function [Gmap, psi] = SENSE_simulation
%
% (c) 2009 - Luis Hernandez-Garcia @ UM
% Last Edit Dec. 15, 2009
%
% do a simulated SENSE reconstruction on cartesian data
close all

N=32;
Npix = N^2;
Ncoils = 8;
x = phantom(N);
x = x(:);

%allS = read_BFields('metaSims128Mhz/coilnometa',8);
allS = fake_Bfields(N);

%psi = coupling16rungs;
psi = eye(8);

Gmap = makeGmap(allS, psi);
Gmap = reshape(Gmap, N,N);
imagesc(abs(Gmap));
title('g-map from Fully sampled S map');

FF = FourierMat(N);  % 2D FT matrix

% Make Encoding matrix including sensitivity: FFS
% there is one for each coil.
for n=1:Ncoils

    smap = allS(n,:);  % a single s-map
    S = repmat(smap,N^2,1);  % replicate it to match the dimensions

    % multiply times the Fourier encoding matrix. This is the new encoding matrix for each coil
    FFS{n} = FF .* S;
end

% a quick test:
allk = zeros(Ncoils, N^2);
allx2 = zeros(Ncoils, N^2);

% forward problem:   generate the kspace data for each coil:
% then FT it to recon the image
for n=1:Ncoils

    k = FF*x;
    x2 = FF*k;
    x2 = x2(end:-1:1, :);
    subplot(234), imagesc(reshape(x,N,N)); axis square, title('object')
    subplot(235), imagesc(reshape(abs(k),N,N)); axis square,  title('k-space with no sens. ')
    subplot(236), imagesc(reshape(abs(x2),N,N)); axis square, title('no sens. weighted');

    z = FFS{n};
    k = z*x ;
    x2 = FF*k ;
    x2 = x2(end:-1:1, :);

    subplot(231), imagesc(reshape(x,N,N)); axis square, title('object')
    subplot(232), imagesc(reshape(abs(k),N,N)); axis square, title('k-space with sens. weight')
    subplot(233), imagesc(reshape(abs(x2),N,N)); axis square, title('sens. weighted');

    allk(n,:) = k';
    allx2(n,:) = x2';


    pause(0.1)
end


% now do a SENSE recon:
% make the unfolding matrix:

% U = makeUnfoldMat(allS, psi);
%
% x3 = U*allx2';
%
% tmp = sum(x3,1);
% figure
% imagesc(reshape(abs(tmp),32,32));
%
%
% fprintf('\nDone with the unaliased data, press <return> to start the aliasing\n');
% pause

%% Now the aliased version
R = 2;
middle = floor(N/2);
stop1 = ceil(N/(R*2))
start1 = floor(N - N/(R*2) + 1)
N2 = floor(N/R);
N2mid = floor(N2/2);

for n=1:Ncoils

    % downsample kspace, then zero-fill the encoded data (no SENSE)
    k = FF*x;
    k = reshape(k,N,N);

    kdown = zeros(size(k));
    kdown(1:stop1, :) = k(1:R:middle, :);
    kdown(start1:end, :) = k(middle+1:R:end, :);
    kdown = kdown(:);

    x2 = FF*kdown;
    x2 = reshape(x2,N,N);
    x2 = x2/max(x2(:));
    x2 = x2(1:R:end, :);
    x2 = [ zeros((N- N/R)/2, N);
        x2;
        zeros((N- N/R)/2, N);
        ];
    x2 = x2(end:-1:1, :);

    subplot(234), imagesc(reshape(x,N,N)); axis square, title('object')
    subplot(235), imagesc(reshape(abs(kdown),N,N)); axis square, title('k-space with no sens.')
    subplot(236), imagesc(reshape(abs(x2),N,N)); axis square ; title('NO sens. weighted');

    % now the  case with sensitivity encoding
    z = FFS{n};
    k = z*x;
    k = reshape(k,N,N);

    kdown = zeros(size(k));
    kdown(1:stop1, :) = k(1:R:middle, :);
    kdown(start1:end, :) = k(middle+1:R:end, :);
    kdown = kdown(:);

    x2 = FF*kdown;
    x2 = reshape(x2,N,N);
    x2 = x2/max(x2(:));
    x2 = x2(1:R:end, :);
    x2 = [ zeros((N- N/R)/2, N);
        x2;
        zeros((N- N/R)/2, N);
        ];
    x2 = x2(end:-1:1, :);

    subplot(231), imagesc(reshape(x,N,N)); axis square, title('object')
    subplot(232), imagesc(reshape(abs(kdown),N,N)); axis square, title('k-space with sens. weight')
    subplot(233), imagesc(reshape(abs(x2),N,N)); axis square; title('sens. weighted');

    allkdown(n,:) = kdown';
    allx2(n,:) = x2(:);

    pause(0.2)
end

PSI = eye(Ncoils);


% %%

figure
for n=1:Ncoils
    imagesc(reshape(abs(allx2(n,:)),N,N));
    drawnow, pause(0.1)
end

Gmap = ones(1,Npix);
x_sense = zeros(1,Npix);
alias_gap = ceil(N/R);
ytmp2 = zeros(Ncoils, Npix);
warning off
for col=1:N
    for row=(N-N/R)/2 : (N +N/R)/2
        % figure out which pixels contribute to alias and extract those
        % sensitivities.  (aliasing occurs in ydirection)

        p = sub2ind([N, N], row, col);
        row_alias = [];
        for r=-R:R
            row_alias =  [row_alias ; row + r*alias_gap];
        end

        % make sure I'm in the FOV
        row_alias = row_alias(row_alias>0);
        row_alias = row_alias(row_alias<=N);
        col_alias = col*ones(size(row_alias));

        % now get their indices in the array - we note that rows = y-coordinate
        p_alias = sub2ind([N, N], row_alias, col_alias);
        
        % encoding goes like this:  y = S * x
        % need to decode by x = inv(S)*y

        % sensitivities at locations that alias on top of each other
        Stmp = allS(:,p_alias);  
        % the aliasing pattern going forward
        % ytmp(:,p) = Stmp * x(p_alias);
        %
        % try to recover x by inversion of Stmp
        % xtmp = pinv(Stmp)*ytmp(:,p);
        xtmp = pinv(Stmp)*allx2(:,p);
        
        x_sense(p_alias) =  xtmp' ;%+ x_sense(p_alias) ;

        Gmap(p_alias) = makeGmap(allS(:,p_alias), 1);
    end
end
subplot(121)
imagesc(reshape(abs(x_sense), N,N));
subplot(122)
imagesc(reshape((Gmap), N,N));


% %
% x3 = U*allx2';
% % x3 = Udown*allx2';
%
% tmp = sum(x3,1);
% %tmp = diag(x3);
% figure
% imagesc(reshape(abs(tmp),N,N));
% for n=1:Ncoils
%     tmp = allx2(:,n);
%     imagesc(reshape(abs(tmp), N,N));
%     pause(0.5)
% end



% Gmap_down = makeGmap(allS_down, psi);
% Gmap_down = reshape(Gmap_down, (N), (N));
% figure
% imagesc(abs(Gmap_down)); title('g-map from aliased S maps')


