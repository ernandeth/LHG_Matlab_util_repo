function [r1 r2  b1err cbf cbv bat res] = img_mrf_lsq(data_dir)
% function [r1 r2 b1err cbf cbv bat residual] = img_mrf_lsq(data_dir)
%
% this function will fit an MRF time series from tseries_file
% to the model specified in gen_signals_vs_230321
% it will use least squares (lsqnonlin) to find the best r1 and r2
% relaxation parameters
savDir = pwd;
cd(data_dir);

aqparms = read_timing_files_asl3dflex('./')
aqparms.order = 1; % <---- set to 2  JUST TESTING TO SEE IF THE ORDER IS WRONG ---
aqparms.flip= deg2rad(5);

raw = readnii('im_mag.nii');

% raw2 = smoothim(raw,'fwhm',3/64);
% msk = makevarmask(raw2 , 75);
% msk(:) = 1;
% lbview(msk);
% raw = raw.*msk;

%{
raw = raw(...
    10:54, ...
    10:54, ...
    20:48, ...
    :);
%}

sz = size(raw);
msk = makevarmask(raw , 80);
lbview(msk); title('Mask')
%msk(:) = 1;
raw = raw.*msk;

r1 = zeros(size(msk));
r2 = zeros(size(msk));
b1err = zeros(size(msk));

prms.r1 = 0;
prms.r2 = 0;
prms.b1err = 0;

cbf = zeros(size(msk));
cbv = zeros(size(msk));
bat = zeros(size(msk));
res = zeros(size(msk));

if numel(sz)==4
timeseries = reshape(raw, sz(1)*sz(2)*sz(3), sz(4));
elseif numel(sz)==3
    timeseries = reshape(raw, sz(1)*sz(2), sz(3));
end

Npix = size(timeseries,1);
Nframes = size(timeseries,2);
%{
% adjust the size of the params. structure if the N frames is less than
% what's the timing files (ie - the experiment was cut short)
aqparms.t_tags      = aqparms.t_tags(1:Nframes);
aqparms.del1        = aqparms.del1(1:Nframes);
aqparms.del2        = aqparms.del2(1:Nframes);
aqparms.del3        = aqparms.del3(1:Nframes);
aqparms.labelcontrol= aqparms.labelcontrol(1:Nframes);
aqparms.doArtSup    = aqparms.doArtSup(1:Nframes);
aqparms.RO_time     = aqparms.RO_time(1:Nframes);
aqparms.t_aq        = aqparms.t_aq(1:Nframes);
%aqparms.RO_time     = aqparms.RO_time*ones(Nframes,1);
%aqparms.t_aq        = aqparms.t_aq*ones(Nframes,1);
%aqparms.order       = 1;
%}

fprintf('\nThese are the acquisition parameters. '); % Begin estimation? (hit ENTER or CTL-C)')
aqparms

aqparms.Ma = [];

parms.f= 0.02;
parms.Mtis0 = 1;
parms. cbva = 0.03;
parms.bat = 0.1;
parms.r1tis= 1;
parms.r2tis = 0.1;
parms.b1err=0;

% generate test data and get the arterial input function
doSub = 0;
dofigs = 1;
figure(33)
%data = gen_signals_vs_230718(parms, aq_parms, dofigs,doSub);
[data Mart]= gen_signals_vs_230918(parms, aqparms, dofigs,doSub);
aqparms.Ma = [];

% Let's recycle the relaxation estimates 
skip_rlx=0;
if skip_rlx
 load estimates.mat 
end

parfor p =1:Npix
    if msk(p) == 1

        % step 1: estimate relaxation parms
        if skip_rlx
            % If we already estimated the relaxation parms, we use them
            a = r1(p);
            b = r2(p);
            c = b1err(p);
        else
            [a b c] = mrf_lsq_rlx(timeseries(p,:), aqparms, [0.01 0.02 0.5], 0);
            showMe = 0;
            r1(p) = a;
            r2(p) = b;
            b1err(p) = c;
        end

        %step 2: estimate vascular parms
        %
        showFit=0;
        if mod(p,1000)==0  % for debugging purposes, show the fit every XX voxels
            fprintf('\nVoxel %d ... %0.2f percent of voxels done',p, 100*p/Npix );
            showFit=1;
            clf
        end
        
        [est R] = mrf_lsq_flow(timeseries(p,:), aqparms, [a b c], showFit);
       
        cbf(p) = est(1);
        cbv(p) = est(2);
        bat(p) = est(3);
        res(p) = R;

        if showFit==1
            title(sprintf('\nVoxel %d ... %0.2f percent of voxels done',p, 100*p/Npix ));
            drawnow
        end

        %}


    end
end

r1 = reshape(r1,size(msk));
r2 = reshape(r2,size(msk));
b1err = reshape(b1err,size(msk));

r1(r1==0)=nan;
r2(r2==0)=nan;


cbf = reshape(cbf,size(msk));
cbv = reshape(cbv,size(msk));

cbf(cbf==0)=nan;
cbv(cbv==0)=nan;

save estimates.mat  r1 r2 b1err cbf cbv res bat
%%
%
doResMask = 0;
resmask = 1;
if doResMask
    resmask = zeros(size(res));
    th = 0.01*std(res(:))
    resmask(res < th) = 1;
end
%figure(66)
% lbview(resmask)

figure(1)
lbview(resmask.*1./r1), title('T1')
lbview(1./r1), title('T1')
caxis([0 2])
colormap parula

figure(2)
lbview(resmask.*1./r2), title('T2')
lbview(1./r2), title('T2')
caxis([0.0 0.15])
colormap parula

figure(3)
lbview(resmask.*b1err), title('B1 error')
lbview(b1err), title('B1 error')
caxis([-1 1]*0.25)
colormap jet

figure(4)
lbview(resmask.*cbf), title('CBF')
lbview(cbf), title('CBF')
caxis([0 0.001])
colormap hot

figure(5)
lbview(resmask.*cbv), title('CBV_a')
lbview(cbv), title('CBV')
caxis([0.0 0.03])
colormap hot

figure(6)
lbview(resmask.*bat), title('BAT')
lbview(bat), title('BAT')
caxis([0.20 0.3])
colormap hot

figure(7)
lbview(resmask.*res), title('Residuals')
lbview(res), title('Residuals: ', sum(res(:)))
caxis([0 1]*1e-2)
colormap gray
%%
%}
return