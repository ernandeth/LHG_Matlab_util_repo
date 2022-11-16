function oCRLB = optCRLB_rand7(params,timing_parms,Nframes,doSub)
Nparms = 7; % 4 for single 5 for two compartment
del = 0.025;
global nobat_tis
sigLen = ( (doSub*Nframes/2) + ((1-doSub)*Nframes) );
D = zeros(Nparms,sigLen);
tags = timing_parms.t_tag;

for par = 1:Nparms
    D(par,:)=differ_opt_rand_v2_7(params,del,par,tags,timing_parms,doSub); % 'differ' for crlb, 'differ_norm' for norm test, 'differ_opt_single' for opt single
end

if nobat_tis
    one_bat = boolean([1 0 1 1 1 1 1]); % no bat_tis in this case
    D = D(one_bat',:);
end
I_theta = D*D';

%     %%%%%%%%%%%
%     std = sqrt(diag(I_theta));
%     N_I_theta = (1./(std*std')).*I_theta;
%
%     figure
%     imagesc(N_I_theta); axis image; caxis([-1 1]); colorbar
%     title('Normalized Fisher Info Matrix')
%     %%%%%%%%%%%

oCRLB = sqrt(diag(inv(I_theta)/10000));
end