function flex_search_150207(dict, parms, raw, msk, searchparms, t1map, betamap, timing_parms)
% function flex_search_150207(dict, parms, raw, msk, searchparms, t1map, betamap, timing_parms)
%
% this function searches over a dictionary
% it uses pixel mask to save time (msk)
%
%  like this:
%           flex_search_150207(dict, parms, raw, msk, [], [], [], [])
% OR ...
%
% if we already know t1map and flip angle, then we use it and generate a
% new dictionary on the fly for each pixel using searchparms
% and the contents of t1map and beta
%   like this:
%           flex_search_150207(dict, parms, raw, msk,searchparms, t1map, betamap, timing_parms )
%

%load dictionary.mat
dict = dict';

% allocate space for resulting maps
xdim = sqrt(size(raw,2));

kfor = zeros(xdim, xdim);
flows = zeros(xdim, xdim);
vols = zeros(xdim, xdim);
trans = zeros(xdim, xdim);
disps = zeros(xdim, xdim);
Rmap = zeros(xdim, xdim);
outbeta = zeros(xdim, xdim);
outt1map = zeros(xdim, xdim);

for p=1: size(raw, 2)
    
    fprintf('\rProgress: %d   out of   %d   pixels', p, xdim*xdim);
    
    if msk(p)==1
        
        % if we already know T1 and flip angle, then we use them.
        % and generate a new dictionary for this pixel
        if ~isempty(searchparms)
            
            searchparms.r1tis = 1/t1map(p);
            searchparms.beta = deg2rad(betamap(p));
            
            [dict, parms] = gen_flex_dictionary_150207 (timing_parms, searchparms);
            dict = dict';
        end
        
        y = raw(:,p);
        %y = detrend(y);
        y = y - mean(y);
        y = (y)/norm(y);
        
        if ~sum(isnan(y))
            [best, R] = best_match(dict, y);
            bp = parms(best);
        else
            bp = parms(1);
            R = 0;
            best = 1;
        end;
        
        if ~isempty(best)
        
            if length(best) >= 2
                fprintf('\rWhoa ... there are %d best matches  ....   !', length(best));
            end
            
            flows(p) = bp(1).f;
            vols(p) = bp(1).cbva;
            trans(p) = bp(1).transit;
            disps(p) = bp(1).Disp;
            outt1map(p) = 1/bp(1).r1tis;
            Rmap(p) = R(1);
            outbeta(p) = rad2deg(bp(1).beta);
            kfor(p) = bp(1).kfor;
        end
    end
    
end


flows =  reshape(flows, xdim, xdim, 1);
trans = reshape(trans, xdim, xdim, 1);
disps = reshape(disps, xdim, xdim, 1);
outt1map = reshape(outt1map, xdim, xdim, 1);
Rmap = reshape(Rmap, xdim, xdim, 1);
vols = reshape(vols, xdim, xdim, 1);
outbeta = reshape(outbeta, xdim, xdim, 1);
kfor = reshape(kfor, xdim, xdim, 1);

beta = outbeta;
t1map = outt1map;

save flows.mat flows
save trans.mat trans
save disps.mat disps
save R_flows.mat Rmap
save kfors.mat kfor
save beta.mat beta
save t1map.mat t1map
save cbva.mat vols

%%
subplot(321)
lightbox(flows * 6000, [],1);
title('Perfusion in ml/100g/min')

subplot(322)
lightbox(trans);
title('transit time (s)')

%
subplot(323)
lightbox((beta), [ ], 1);
title('Flip Angle (deg)')

subplot(324)
lightbox(Rmap);
title('R map');

subplot(325)
lightbox(t1map);
title('T1 map  (s)');

subplot(326)
lightbox(Rmap, [], 1);
title('Best Correlation (R)');

subplot(326)
lightbox(vols, [], 1);
title('CBV fraction');


return
