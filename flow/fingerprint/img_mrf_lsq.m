function [r1 r2  b1err cbf cbv bat res] = img_mrf_lsq(tseries_file)
% function [r1 r2 b1err cbf cbv bat residual] = img_mrf_lsq(tseries_file)
%
% this function will fit an MRF time series from tseries_file
% to the model specified in gen_signals_vs_230321
% it will use least squares (lsqnonlin) to find the best r1 and r2
% relaxation parameters

aqparms = read_timing_files('./')
%aqparms.order = 2; % <----  JUST TESTING TO SEE IF THE ORDER IS WRONG ---
if(aqparms.RO_type=='GRE')
    prms.flip = deg2rad(20);
end

raw = readnii(tseries_file);

raw2 = smoothim(raw,'fwhm',3/64);
msk = makevarmask(raw2 , 75);
msk(:) = 1;
lbview(msk); 
raw = raw.*msk;

%{
raw = raw(...
    10:54, ...
    10:54, ...
    20:48, ...
    :);
%}

sz = size(raw);
%msk = makevarmask(raw , 55);
%lbview(msk)
%msk(:) = 1;
%raw = raw.*msk;

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

timeseries = reshape(raw, sz(1)*sz(2)*sz(3), sz(4));
Npix = size(timeseries,1);
Nframes = size(timeseries,2);

% adjust the size of the params. structure if the N frames is less than
% what's the timing files (ie - the speriment was cut short)
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

aqparms.order       = 1;
fprintf('\nThese are the acquisition parmeteres.  Begin estimation? (hit ENTER or CTL-C)')
aqparms
pause

parfor p =1:Npix
    if msk(p) == 1
        fprintf('\nVoxel %d ... %0.2f percent of voxels done',p, 100*p/Npix );

        % step 1: estimate relaxation parms
        [a b c] = mrf_lsq_rlx(timeseries(p,:), aqparms, 0);
        showMe = 0;
        r1(p) = a;
        r2(p) = b;
        b1err(p) = c;

        %step 2: estimate vascular parms
        %
        [est R] = mrf_lsq_flow(timeseries(p,:), aqparms, [a b c], 0);
        cbf(p) = est(1);
        cbv(p) = est(2);
        bat(p) = est(3);
        res(p) = R;
        %}
            

        if mod(p,500)==0  % for debugging purposes, show the fit every 1000 voxels
            [est R] = mrf_lsq_flow(timeseries(p,:), aqparms, [a b c], 1);
        end

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
doResMask = 1;
resmask = 1;
if doResMask
    resmask = zeros(size(res));
    th = 0.05*std(res(:))
    resmask(res < th) = 1;
end
lbview(resmask)
%%
figure(1)
orthoview(resmask.*1./r1), title('T1')
lbview(1./r1), title('T1')
caxis([0 3.5])
colormap parula

figure(2)
orthoview(resmask.*1./r2), title('T2')
lbview(1./r2), title('T2')
caxis([0.0 0.15])
colormap parula

figure(3)
orthoview(resmask.*b1err), title('B1 error')
lbview(b1err), title('B1 error')
caxis([-1 1]*0.05)
colormap jet

figure(4)
orthoview(resmask.*cbf), title('CBF')
lbview(cbf), title('CBF')
caxis([0 0.015])
colormap hot

figure(5)
orthoview(resmask.*cbv), title('CBV_a')
lbview(cbv), title('CBV')
caxis([0.0 0.025])
colormap hot

figure(6)
orthoview(resmask.*bat), title('BAT')
lbview(bat), title('BAT')
caxis([0.0 0.5])
colormap hot

figure(7)
orthoview(resmask.*res), title('Residuals')
lbview(res), title('Residuals')
caxis([0 1]*1e-2)
colormap hot
%%
%}
return