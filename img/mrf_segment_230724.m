function [tcomponents class] = mrf_segment_230724(filename, Ncomps, NITER)
%function tcomponents = mrf_segment(filename, Ncomps or T1 values....)
%

if nargin<2
    Ncomps = 3;
end

initSVD = 0   % it seems to be better not to initialize with PCA.  Maybe physio components dominate

r1_array = linspace(5,0.2,Ncomps);
%r2_array = linspace(3,20,Ncomps);
r2_array = linspace(15,15,Ncomps);  % make it only T1 dependent
r1_array = 1./linspace(0.8,3,Ncomps);

%r2_array = 1./[0.04 0.09 0.25] ;  % make it only T1 dependent
%r1_array = 1./[0.9  1.4   3];

[raw hdr] = readnii(filename);
% reshape data into a carpet plot
xdim = size(raw,1);
ydim = size(raw,2);
zdim = size(raw,3);
tdim = size(raw,4);

raw = reshape(raw, [xdim*ydim*zdim , tdim] )';

% make a mask based on the variance
TopPercentSigma = 55;
sigma = std(raw,[],1);
msk = ccvarmask(sigma, TopPercentSigma);
%{
figure;
subplot(211)
lbview(reshape(sigma, [xdim, ydim, zdim]));
title('std. dev. map')
subplot(212)
lbview(reshape(msk, [xdim, ydim, zdim]));
title('mask for analysis');
%}
mraw = raw .* repmat(msk, tdim,1);
mraw(isinf(mraw))=0;
mraw(isnan(mraw))=0;

Npix = size(mraw,2);
innerprod = zeros(1,Ncomps);

% initialise base signals:
tcomp_init = rand(tdim, Ncomps);

parms.f         = 50/6000;
parms.Mtis0     = 1;
parms. cbva     = 0.02;
parms.bat       = 0.1;
parms.r1tis     = 1/1.4;
parms.r2tis     = 1/0.05;
parms.flip      = deg2rad(90);
parms.b1err     = 0;

aq_parms = read_timing_files('./')

% initializing components' timcourses.
% calculate what the original components are
if initSVD
    [u,s,v] = svd(mraw,"econ");
end
%
for n=1:Ncomps
    if initSVD
        tcomp_init(:,n) = u(:,n);
    else
        parms.r1tis = r1_array(n);
        parms.r2tis =  r2_array(n);
        obs = gen_signals_vs_230718(parms, aq_parms, 0,0);
        obs = obs(1:tdim);
        tcomp_init(:,n) = obs';
    end
    % alternatively, we can start with random conponents (does not do well)
    %tcomp_init = randn(size(tcomp_init));
end

%}

%{
subplot(211)
plot(tcomp_init)
title( 'seed components from simulation')
drawnow
%}

%{
inds = find(msk);
tmp = mraw(:,inds);
[idx c] = kmeans(tmp',Ncomps);
tcomp_init= c';
%}

% begin classification.
% The choices start out as the initial possible components
% They will evolve as we iterate and add timecourses.

%tcomp_init = orth(tcomp_init);

for n=1:Ncomps
    tcomp_init(:,n) = tcomp_init(:,n) - mean(tcomp_init(:,n));
    tcomp_init(:,n) = tcomp_init(:,n) / norm(tcomp_init(:,n));
end
tcomponents = tcomp_init;


tcom_init = tcomp_init(20:100,:);
tcomponents=tcomponents(20:100,:);
mraw = mraw(20:100,:);
raw = raw(20:100,:);

tdim=81;

for m=1:NITER

    class = zeros(size(msk));
    for p=1:Npix
        maxnorm = 0;

        if msk(p)==1;

            % zero mean and normalize each timecourse
            tmp2 = mraw(:,p) ;
            tmp2 = tmp2-mean(tmp2);
            tmp2 = tmp2/ norm(tmp2);

            % find which seed component is most correlated with this time course.
            for c=1:Ncomps
                tmp = tcomponents(:,c); 
                innerprod(c) = tmp' * tmp2;
                %cc=corrcoef(tmp,tmp2);
                %innerprod(c) = cc(1,2);
            end

            [junk cmax] = max((innerprod));

            % assign this pixel to a class (corresponds to one of seed timecourse)
            class(p) = cmax;

            % update the components:
            % add timecourse to best candidate component:
            tcomponents(:, cmax) =   tcomponents(:, cmax) +  tmp2;
            % zero mean and normalize the new component
            %tcomponents(:, cmax) =   tcomponents(:, cmax) - mean(tcomponents(:, cmax));
            tcomponents(:, cmax) =   tcomponents(:, cmax) / norm(tcomponents(:,cmax));
            %plot(tcomponents); drawnow
        end
    end


end


% sort the components according to the amount of energy
nrg = sum(tcomponents.^2, 1);
[nrg2 inds] = sort(nrg,'descend');
tmp = tcomponents(:,inds);
tcomponents = tmp;

% force the components to span from 0 to 1
for n=1:Ncomps
    tcomponents(:, n) =   tcomponents(:, n) - min(tcomponents(:, n));
    tcomponents(:, n) =   tcomponents(:, n) / max(tcomponents(:,n));
    %plot(tcomponents); drawnow
end
% Make a model of the signal (design matrix)
q = ([ones(tdim,1) tcomponents ]);
% perform regression of model on the data
fractions = pinv(q)*raw;

% we don't care about the DC baseline
fractions = fractions(2:end,:);

% calculate the amount of signal accounted for each component
Total = sum(fractions , 1);
fractions = fractions ./ Total;
fractions = fractions .* msk;
fractions = reshape(fractions', [xdim ydim zdim Ncomps]);
writenii('beta_fractions', fractions);

figure
crop = round(xdim/16);
tmp = fractions(crop+1:end-crop, crop+1:end-crop,crop+1:end-crop, :);
for n=1:Ncomps
    subplot(2, ceil(Ncomps/2),  n)
    lbview(tmp, 'frame',n); %colorbar off
    title('fractional content')
    caxis([-1 1])
end
colormap jet

% show the classification
class = reshape(class, [xdim ydim zdim]);
tmp = class(crop+1:end-crop, crop+1:end-crop,crop+1:end-crop, :);
figure
lbview(tmp); title('Voxel Classifications')

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
