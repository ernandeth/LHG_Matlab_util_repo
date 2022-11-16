function flex_search_150208(dict, parms, raw, msk, priors,  timing_parms)
% function flex_search_150208(dict, parms, raw, msk, priors, timing_parms)
%
% this function searches over a dictionary.
% inputs
%       dict = dictionary to use in the search
%       parms = corresponding parms to each dictionary entry
%       raw = input data (image time series as a matrix : frames x pixels
%       msk = binary mask so you can skip voxels without content
%
%  like this:
%           flex_search_150208(dict, parms, raw, msk, [], timing_parms)
%
% OR ...
%
% if we already know some parm maps, like t1map and flip angle, then we use them and generate a
% new dictionary on the fly for each pixel using the priors structure.
%   like this:
%
%           flex_search_150208([], [], raw, msk, priors, timing_parms )
%
%  - You can tell priors a range for the dictionary or
%  - You let the program know a known parameter maps by sticking in a file name for
%       the .mat file that contains them.
%  - it ignores the dict and parms in this case
%
% For example:
%
% priors.kfor = 0.02;
% priors.m0tis = 1
% priors.r1tis = 't1map.mat'    % <-----
% priors.beta = 'flipang.mat'    % <-----
% priors.f = linspace(0,100,20)/6000;
% priors.cbva =      0.01 ;
% priors.transit=   linspace(0.5, 4, 10);
% priors.L = 1;
% priors.Disp =      20;
% priors.Ptime =     0.5;
%
% Note that you also have to tell the prog. the timing parms for the acquisition so that it
% knows how to make the dictionaries:
%
% For example:
% timing_parms.PID = ones(10,1)*1.2;
% timing_parms.label_duration = ones(10,1)*1.6;
% timing_parms.t_adjust = ones(10,1);
%


%load dictionary.mat
dict = dict';

% allocate space for resulting maps
xdim = sqrt(size(raw,2));

kformap = zeros(xdim, xdim);
flowmap = zeros(xdim, xdim);
volmap = zeros(xdim, xdim);
transmap = zeros(xdim, xdim);
dispmap = zeros(xdim, xdim);
Rmap = zeros(xdim, xdim);
flipmap = zeros(xdim, xdim);
t1map = zeros(xdim, xdim);

% load the knowns maps from file
if ~isempty(priors)
    if ischar(priors.f)
        load(priors.f);
        flowmap = flowmap(:);
    end
    if ischar(priors.r1tis)
        load(priors.r1tis);
        t1map = t1map(:);
    end
    if ischar(priors.transit)
        load(priors.transit);
        transmap = transmap(:);
    end
    if ischar(priors.beta)
        load(priors.beta);
        flipmap = flipmap(:);
    end
    if ischar(priors.cbva)
        load(priors.cbva);
        volmap = volmap(:);
    end
    if ischar(priors.kfor)
        load(priors.kfor);
        kformap = kformap(:);
    end
    dparms = priors;
end

Npix = size(raw,2);

for p=1: size(raw, 2)
    
    %fprintf('\rProgress: %d   out of   %d   pixels', p, xdim*xdim);
    
    if msk(p)==1 || size(raw,2)==1
        
        if ~isempty(priors)
            
            % if we already know T1 and flip angle, then we use them.
            % and generate a new dictionary for this pixel.
            % otherwise, we use the same dictionary for all pixels.
            if ischar(priors.f)
                dparms.f = flowmap(p);
            end
            if ischar(priors.r1tis)
                dparms.r1tis = 1/t1map(p);
            end
            if ischar(priors.transit)
                dparms.transit = transmap(p);
            end
            if ischar(priors.beta)
                dparms.beta = deg2rad(flipmap(p));
            end
            if ischar(priors.cbva)
                dparms.cbva = volmap(p);
            end
            if ischar(priors.kfor)
                dparms.kfor = kformap(p);
            end
            
            [dict, parms] = gen_flex_dictionary_150207 (timing_parms, dparms);
            dict = dict';
            
        end
        
        
        y = raw(:,p);
        %y = detrend(y);
        y = y - mean(y);
        y = (y)/norm(y);
        
        %{
        % this is intended to normalize and mean center  each half of the
        % data separately for the cases where we're concatenating two runs
        
        y = normalize_halves(y);
       
       for n=1:size(dict,2)
           dict(:,n) = normalize_halves (dict(:,n)) ;
       end
       %}    
        
        
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
            
            flowmap(p) = bp(1).f;
            volmap(p) = bp(1).cbva;
            transmap(p) = bp(1).transit;
            dispmap(p) = bp(1).Disp;
            t1map(p) = 1/bp(1).r1tis;
            Rmap(p) = R(1);
            flipmap(p) = rad2deg(bp(1).beta);
            kformap(p) = bp(1).kfor;
        end
        
        
        if Npix==1
            fprintf('\nROI analysis .... \n');
            close all;
            figure
            plot(dict);
            hold on
            plot(y,'k')
            plot(y,'k*')
            fprintf('\nEstimated parms:   \n');
            bp(1)
            
            drawnow
            pause(1)
        end
        
    end
    
end


Rmap(isnan(Rmap)) = 0;

if Npix>1
    flowmap =  reshape(flowmap, xdim, xdim, 1);
    transmap = reshape(transmap, xdim, xdim, 1);
    dispmap = reshape(dispmap, xdim, xdim, 1);
    t1map = reshape(t1map, xdim, xdim, 1);
    Rmap = reshape(Rmap, xdim, xdim, 1);
    volmap = reshape(volmap, xdim, xdim, 1);
    flipmap = reshape(flipmap, xdim, xdim, 1);
    kformap = reshape(kformap, xdim, xdim, 1);
end

save flows.mat flowmap
save trans.mat transmap
save disps.mat dispmap
save R_flows.mat Rmap
save kfors.mat kformap
save flipang.mat flipmap
save t1map.mat t1map
save cbva.mat volmap

%%
if Npix>1
    figure
    subplot(321)
    lightbox(flowmap * 6000, [],1);
    title('Perfusion in ml/100g/min')
    
    subplot(322)
    lightbox(transmap);
    title('BAT (s)')
    
    %
    subplot(323)
    lightbox((flipmap), [ ], 1);
    title('Flip Angle (deg)')
    
    subplot(324)
    lightbox(kformap);
    title('MT rate (1/s)');
    
    subplot(325)
    lightbox(t1map);
    title('T1 map  (s)');
    
    subplot(326)
    lightbox(volmap, [0 0.03], 1);
    title('CBV fraction');
    
    figure
    lightbox(Rmap,[0.5 1],1);
    title('Best match score (R)')
end

return

function result= normalize_halves(y)

 %% kludge for concatenated time courses:
        y1 = y(1:end/2);
        y2 = y(end/2+1:end);

        y1 = y1 - mean(y1);
        y1 = (y1)/norm(y1);
        
        y2 = y2 - mean(y2);
        y2 = (y2)/norm(y2);
        
        result = [y1 ; y2];
 return
        
