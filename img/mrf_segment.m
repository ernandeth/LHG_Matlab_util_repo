function tcomponents = mrf_segment(filename, Ncomps, NITER)
%function tcomponents = mrf_segment(filename, Ncomps or T1 values....)
%
if length(Ncomps)>1;
    r1_array = 1./Ncomps
    Ncomps = length(Ncomps);
else
    r1_array = linspace(3,0.3,Ncomps);
end

if nargin<2
    Ncomps = 3;
end


[raw hdr] = read_img(filename);
if isfield(hdr,'originator')
    hdr=nii2avw_hdr(hdr);
end

TopPercentSigma = 80;
sigma = std(raw,1);
msk = ccvarmask(sigma, TopPercentSigma);
%{
figure;
subplot(211)
lightbox(reshape(sigma,hdr.xdim, hdr.ydim, hdr.zdim));
title('std. dev. map')
subplot(212)
lightbox(reshape(msk,hdr.xdim, hdr.ydim, hdr.zdim));
title('mask for analysis');
%}
mraw = raw .* repmat(msk, hdr.tdim,1);
mraw(isinf(mraw))=0;
mraw(isnan(mraw))=0;

Npix = size(mraw,2);
innerprod = zeros(1,Ncomps);

% initialise base signals:
seeds = rand(hdr.tdim, Ncomps);

parms.f         = 50/6000;
parms.Mtis0     = 1;
parms. cbva     = 0.02;
parms.bat       = 0.35;
parms.r1tis     = 1/1.4;
parms.flip      = deg2rad(90);
parms.alpha_ai  = 0.8; % arterial inversion efficiency
parms.alpha_ti  = 0.75; % tissue inversion efficiency
parms.alpha_ts  = 0.17; % T2 weighting in tissue due to arterial suppression

tmp = load('t_delays.txt');
aq_parms.t_delays       = load('t_delays.txt');
aq_parms.t_adjusts      = load('t_adjusts.txt');
aq_parms.labelcontrol   = load('isVelocitySelective.txt');
aq_parms.doArtSup       = load('doArtSuppression.txt');
aq_parms.ArtSup_delay   = 0.1 *ones(hdr.tdim,1); % delay between AS pulse and acqusition
aq_parms.t_aq           = 0.58* ones(hdr.tdim,1);
aq_parms.t_tag          = zeros(hdr.tdim,1);

aq_parms.t_tag(:) = 0;
% some of the acquisitions got clipped!)
aq_parms.t_delays = aq_parms.t_delays(1:hdr.tdim);
aq_parms.t_adjusts = aq_parms.t_adjusts(1:hdr.tdim);
aq_parms.labelcontrol = aq_parms.labelcontrol(1:hdr.tdim);
aq_parms.doArtSup = aq_parms.doArtSup(1:hdr.tdim);

%
for n=1:Ncomps
    parms.r1tis = r1_array(n);
    obs = gen_signals_vs_201020(parms, aq_parms, 0,0);
    seeds(:,n) = abs(obs'/obs(1));
end

tcomponents = seeds;

%{
subplot(211)
plot(seeds)
title( 'seed components from simulation')
drawnow
%}

%{
[idx c] = kmeans(mraw',Ncomps);
tcomponents = [c' ones(hdr.tdim,1)];
%}

class = zeros(size(msk));
for m=1:NITER
    for p=1:Npix
        maxnorm = 0;
        %mraw(:,p) = mraw(:,p)/mraw(1,p); % scale
        %mraw(:,p) = mraw(:,p)/norm(mraw(1,p)); % scale
        
        if msk(p)==1;
            % find which component is most correlated with this time course.
            for c=1:Ncomps
                tmp = tcomponents(:, c) / norm(tcomponents(:,c));
                innerprod(c) = tmp' * mraw(:,p);
            end
            
            [junk cmax] = max(innerprod);
            %[innerprod cmax];
            class(p) = cmax;
            % add timecourse to best candidate:
            tcomponents(:, cmax) =   tcomponents(:, cmax) +  mraw(:,p);
        end
    end
end



% pick the top 3 components (with most energy)
nrg = zeros(3,1);
for n=1:Ncomps
    nrg(n) = norm(tcomponents(:,n),2);
end
[nrg2 inds] = sort(nrg,'descend');
nrg2
tcomponents = tcomponents(:,inds(1:Ncomps));
%Ncomps = 3;
%
% reorder them so that they approximate the order of the initial seeds:
% figure
order = zeros(Ncomps,1);
for n=1:Ncomps
    tcomponents(:,n) = abs(tcomponents(:,n)/ norm(tcomponents(:,n)));
    seeds(:,n) = abs(seeds(:,n)/ norm(seeds(:,n)));
    
    % check for the largest correlation with the seeds
    R = zeros(Ncomps,1);
    for m=1:Ncomps
        c = corrcoef(seeds(50:end,m),  tcomponents(50:end,n))  ;
        R(m) = c(1,2);
        R(m) = tcomponents(end,m);
    end
%     subplot(Ncomps,1,n)
%     plot(seeds); hold on; legend('1')
%     plot(tcomponents(:,n), 'k')
    order(n) = find(R == max(R));
end
[junk order] = sort(tcomponents(end,:));

sprintf('\nOrder of the components: %d', order)
tmp = zeros(size(tcomponents));
for n=1:Ncomps
    tmp(:,n) = tcomponents(:,order(n));
end
tcomponents = tmp;
% 
tcomponents = [tcomponents ones(hdr.tdim,1)];
contrast = eye(Ncomps+1);
contrast(:, end) = -1;

spmJr(filename, abs(tcomponents), contrast);
%{
for n=1:Ncomps+1
    figure
    subplot(211)
    plot(tcomponents(:,n)), title (sprintf('component %d',n));
    set(gcf,'Name', 'time components from clustering')
    hold on
    
    subplot(212)
    str = sprintf('Zmap_%04d.img', n);
    z = lightbox(str);
    lightbox(z .* reshape(msk,size(z)));
end
%}
[fractions h] = read_img('Bhats.img');
fractions = fractions(1:end-1,:); 
Total = sum(fractions , 1);
fractions = fractions ./ repmat(Total, Ncomps,1);
fractions = fractions .* repmat(msk, Ncomps,1);
h.tdim = h.tdim-1;
write_img('beta_fractions.img', 1e4*fractions, h);

figure
crop = round(hdr.xdim/12);

for n=1:Ncomps
    subplot(1, Ncomps,  n)
    tmp = lightbox('beta_fractions.img',[0 1e4],[],n);
    tmp = flip(tmp(crop+1:end-crop, crop+1:end-crop,:),2);
    lightbox(tmp/100, [0 1e2],[]); %colorbar off
end
colormap parula
%{
figure
for n=1:Ncomps-1
    subplot(Ncomps-1,1,n); plot(tcomponents(:,n)/tcomponents(1,n)); hold on; plot(seeds(:,n)/seeds(1,n));
    legend('final', 'initial')
end
%}
return

function msk = ccvarmask(varimage , th)
% makes a mask that preserves the data with the top TH percentage of the
% values

ordered = sort(varimage(:));
Nintgrl = cumsum(ordered)/sum(ordered(:)) * 100;
thval = find(Nintgrl>100-th);
%subplot(211), plot(ordered); title('Ordered Std. Deviations')
%subplot(212), plot(Nintgrl); title('Intrageted , Normalized Std. Deviations')

thval = ordered(thval(1))
msk = varimage;
msk(msk < thval) = 0;
msk(msk>=thval) = 1;
return
