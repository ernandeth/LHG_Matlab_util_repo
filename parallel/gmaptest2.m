% testing sense/g-maps

Np = 64;
R = 3;
Nc = 8;

alspace = ceil(Np/R);
extra = alspace*R-Np;
arp = [-Np/2:Np/2-1+extra];
Gmap = ones([Np alspace*R]);

%ark = [-Np/2:R:Np/2-1];
%Nk = length(ark);

psid = repmat(1:Nc,Nc,1) - repmat([1:Nc]',1,Nc);
psi = (.1).^abs(psid);
%Gaussian Smaps
sctr = linspace(-Np/2,Np/2-1,Nc);
arpa = repmat(arp,Nc,1);
sctra = repmat(sctr',1,Np+extra);
Smaps = exp(-(arpa - sctra).^2/50);

% Folding array
F = ones(size(exp(i*2*pi*[0:R-1]/R)));
for ploop = 1:alspace
    SE = Smaps(:,ploop:alspace:alspace*R).*repmat(F,Nc,1);
    % Gmap 1D
    part1 = inv(SE'*SE);
    part2 = SE'*SE;
    %part2 = Smaps(smrep,:)'*Smaps(smrep,:);
    Gmaptmp = real(sqrt(diag(part1).*diag(part2)));
    Gmap1(ploop:alspace:alspace*R) = Gmaptmp;
end


%OK now in 2D
for xloop = 1:Np
    for ploop = 1:alspace
        yloc = ploop:alspace:alspace*R;
        keep = (sqrt((yloc-(Np/2+1)).^2 + (xloop-(Np/2+1)).^2) <= Np/2+1);
        % keep is a circle of Np/2 radius - this is where sensitvity data
        % are
        
        SE = Smaps(:,yloc(keep)).*repmat(F(keep),Nc,1);
        % Gmap 1D
        part1 = inv(SE'*SE);
        part2 = SE'*SE;
        %part2 = Smaps(smrep,:)'*Smaps(smrep,:);
        Gmaptmp = real(sqrt(diag(part1).*diag(part2)));
        Gmap(xloop,yloc(keep)) = Gmaptmp;
    end
end
imagesc(Gmap);colormap(gray);axis('image'); colorbar
