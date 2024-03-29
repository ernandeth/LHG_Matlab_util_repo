function [tcomponents class] = mrf_segment(filename, Ncomps, NITER)
%function tcomponents = mrf_segment(filename, Ncomps or T1 values....)
%
if length(Ncomps)>1;
    r1_array = 1./Ncomps
    Ncomps = length(Ncomps);
else
    r1_array = linspace(3,0.25,Ncomps);
    r2_array = linspace(60,3,Ncomps);
end

if nargin<2
    Ncomps = 3;
end


[raw hdr] = read_img(filename);
if isfield(hdr,'originator')
    hdr=nii2avw_hdr(hdr);
end

TopPercentSigma = 50;
sigma = std(raw,1);
msk = ccvarmask(sigma, TopPercentSigma);
%
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
tcomp_init = rand(hdr.tdim, Ncomps);

parms.f         = 50/6000;
parms.Mtis0     = 1;
parms. cbva     = 0.02;
parms.bat       = 0.1;
parms.r1tis     = 1/1.4;
parms.r2tis     = 1/0.05;
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
aq_parms.t_delay = aq_parms.t_delays(1:hdr.tdim);
aq_parms.t_adjusts = aq_parms.t_adjusts(1:hdr.tdim);
aq_parms.labelcontrol = aq_parms.labelcontrol(1:hdr.tdim);
aq_parms.doArtSup = aq_parms.doArtSup(1:hdr.tdim);

aq_parms.readout_type = 'FSE';
aq_parms.label_type = 'BIR8inv';
aq_parms.order=1;
%
for n=1:Ncomps
    parms.r1tis = r1_array(n);
    %obs = gen_signals_vs_201020(parms, aq_parms, 0,0);
    parms.r2tis =  r2_array(n);
    obs = gen_signals_vs_230321(parms, aq_parms, 0,0);
    
    % use a random time course to seed the components
    tcomp_init(:,n) = obs';
    tcomp_init(:,n) = tcomp_init(:,n) - mean(tcomp_init(:,n));
    tcomp_init(:,n) = tcomp_init(:,n) / norm(tcomp_init(:,n));

    %tcomp_init = randn(size(tcomp_init));
   
end
tcomponents = tcomp_init;

%{
subplot(211)
plot(tcomp_init)
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
        
        if msk(p)==1;

            % normalize each timecourse                
            tmp2 = mraw(:,p) / norm(mraw(:,p));

            % find which seed component is most correlated with this time course.
            for c=1:Ncomps
                tmp = tcomponents(:, c) / norm(tcomponents(:,c));
                innerprod(c) = tmp' * tmp2;
                %cc=corrcoef(tmp,tmp2);
                %innerprod(c) = cc(1,2);

            end            
            [junk cmax] = max(innerprod);
            if junk>1e-2
                % assign this pixel to a class (corresponds to one of seed timecourse)
                class(p) = cmax;

                % update the components:
                % add timecourse to best candidate seed:
                tcomponents(:, cmax) =   tcomponents(:, cmax) +  tmp2;
                % normalize the seed
                tcomponents(:, cmax) =   tcomponents(:, cmax) /norm(tcomponents(:,cmax));
            end
        end
    end


end


% sort the components according to the amount of energy
nrg = sum(tcomponents.^2, 1);
[nrg2 inds] = sort(nrg,'descend');
tmp = tcomponents(:,inds);
tcomponents = tmp;
tcomponents = tcomponents - mean(tcomponents,1);
tcomponents = tcomponents./abs(max(tcomponents,[],1));
% 
tcomponents = [tcomponents ones(hdr.tdim,1)];
contrast = eye(Ncomps+1);
contrast(end,:) = [];

spmJr(filename, (tcomponents), contrast);
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
[fractions h] = read_img('ConBhats.img');
Total = sum(fractions , 1);
fractions = fractions ./ Total;
fractions = fractions .* msk;
write_img('beta_fractions.img', 1e4*fractions, h);

figure
crop = round(hdr.xdim/12);

for n=1:Ncomps
    subplot(2, ceil(Ncomps/2),  n)
    tmp = lightbox('beta_fractions.img',[0 1e4],[],n);
    tmp = flip(tmp(crop+1:end-crop, crop+1:end-crop,crop+1:end-crop),2);
    lightbox(tmp/100, [-1 1]*100,[]); %colorbar off
end
colormap parula
%
figure
for n=1:Ncomps
    subplot(Ncomps,1,n); 
    plot(tcomp_init(:,n)); hold on;
    plot(tcomponents(:,n)); hold off;
    legend('Initial', 'final')
    hold off
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
