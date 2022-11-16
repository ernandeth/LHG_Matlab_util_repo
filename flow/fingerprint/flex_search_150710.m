function flex_search_150710(dict, parms, raw, msk, priors)
% function flex_search_150710(dict, parms, raw, msk, priors)
%
% this function searches over a dictionary.
% inputs
%       dict = dictionary to use in the search
%       parms = corresponding parms to each dictionary entry
%       raw = input data (image time series as a matrix : frames x pixels
%       msk = binary mask so you can skip voxels without content
%      priors = parms whos values you already know
%
% if we already know some parm maps, like t1map and flip angle, then we use them
% TO MASK OUT ENTRIES IN THE DICTIONARY THAT DO NOT MATCH
% THIS VERSION GENERATES A SINGLE DICTIONARY AND MASKS OUT INELLIGIBLE
% ENTRIES.
% NOTE THAT YOU HAVE TO GIVE IT A COMPLETE DICTIONARY AS AN INPUT
%
%   like this:
%
%           flex_search_150710( DICT, PARMS, raw, msk, priors )
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
% priors.flip = 'flipang.mat'    % <-----
% priors.f = linspace(0,100,20)/6000;
% priors.cbva =      0.01 ;
% priors.bat=   linspace(0.5, 4, 10);
% priors.bat2=  0.5;
% priors.L = 1;
% priors.Disp =      20;
% priors.Ptime =     0.5;
%
% Note that you first need to generate a dictionary
%
%
% For example:
%
% [dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);
% flex_search_150710(dict, parms, obs, [], priors)


%load dictionary.mat
dict = dict';

% allocate space for resulting maps
xdim = sqrt(size(raw,2));

kformap = zeros(xdim, xdim);
flowmap = zeros(xdim, xdim);
volmap = zeros(xdim, xdim);
transmap = zeros(xdim, xdim);
transmap2 = zeros(xdim, xdim);
dispmap = zeros(xdim, xdim);
Rmap = zeros(xdim, xdim);
flipmap = zeros(xdim, xdim);
t1map = zeros(xdim, xdim);

saveflip = 1;
saveflow = 1;
savebat =1 ;
savebat2 =1 ;

savedisps = 1;
savekfors = 1;
savet1 = 1;
savecbv = 1;

% load the knowns maps from file
if ~isempty(priors)
    
    if isfield(priors,'f')   && ischar(priors.f)
        load(priors.f);
        flowmap = flowmap(:);
        %saveflow = 0;
        
    end
    
    if isfield(priors,'r1tis') && ischar(priors.r1tis)
        load(priors.r1tis);
        t1map = t1map(:);
        %savet1 =0;
    end
    
    if isfield(priors,'bat') && ischar(priors.bat)
        load(priors.bat);
        transmap = transmap(:);
        %savebat = 0;
    end
    
    if isfield(priors,'bat2') && ischar(priors.bat2)
        load(priors.bat2);
        transmap2 = transmap2(:);
        %savebat = 0;
    end
    
    if isfield(priors,'flip') && ischar(priors.flip)
        load(priors.flip);
        flipmap = flipmap(:);
        %saveflip = 0;
    end
    
    
    if isfield(priors,'cbva') && ischar(priors.cbva)
        load(priors.cbva);
        volmap = volmap(:);
        %savecbv = 0;
    end
    
    if isfield(priors,'kfor') && ischar(priors.kfor)
        load(priors.kfor);
        kformap = kformap(:);
        %savekfors = 0;
    end
    
    
    dparms = priors;
    
end
%{
% Normalize all the entries of the dictionary
fprintf('\n Normalizing dictionary ...');

for n=1:size(dict,1)
    entry = dict(n,:);
    tmp = entry - mean(entry);
    tmp = tmp / norm(tmp);
    dict(n,:) = tmp;
end
 %}
    
Npix = size(raw,2);
pr = [];
% the dictionary entries must have values within 15% of the value of the prior
matchfactor = 0.15; 
        
parfor p=1: size(raw, 2)
    
    if ~mod(p,200)
        fprintf('\rProgress: %d   out of   %d   pixels', p, Npix);
    end
    
    if msk(p)==1 || size(raw,2)==1
        % searchables is a mask indicating whith entries of the dictionary
        % can be used for this pixel
        searchables = ones(size(dict,2),1);
        
        pr = priors;
        
        % if we already know T1 and flip angle, then we use them.
        % if the priors are maps, then we look at each pixel.
        if isfield(priors,'f')       
            if ischar(priors.f)
                pr.f = flowmap(p);
            end
            
            % get the flow values corresponding to each entry of the
            % dictionary:
            allf = cell2mat({parms(:).f});
            dict_msk = zeros(size(searchables));
            % see which of those flow values in the dictionary are close to
            % the prior known value.  Those entries can be searched.
            inds = find( abs(allf-pr.f) < matchfactor*pr.f);
         
            dict_msk(inds) = 1;
            if ~isempty(inds),
                searchables = searchables .* dict_msk;
            end
            
        end
        
        
        if isfield(priors,'r1tis')
            if ischar(priors.r1tis)
                pr.r1tis = 1/t1map(p);
            end
            
            allr1tis = cell2mat({parms(:).r1tis});
            dict_msk = zeros(size(searchables));
            inds = find( abs(allr1tis-pr.r1tis) < matchfactor*pr.r1tis);
            dict_msk(inds) = 1;
            if ~isempty(inds),
                searchables = searchables .* dict_msk;
            end
        end
        
        if isfield(priors,'bat')
            if ischar(priors.bat)
                pr.bat = transmap(p);
            end
            allbat = cell2mat({parms(:).bat});
            dict_msk = zeros(size(searchables));
            inds = find( abs(allbat-pr.bat) < matchfactor*pr.bat);
            dict_msk(inds) = 1;
            if ~isempty(inds),
                searchables = searchables .* dict_msk;
            end
        end
        
        if isfield(priors,'bat2')
            if ischar(priors.bat2)
                pr.bat2 = transmap2(p);
            end
            allbat2 = cell2mat({parms(:).bat2});
            dict_msk = zeros(size(searchables));
            inds = find( abs(allbat2-pr.bat2) < matchfactor*pr.bat2);
            dict_msk(inds) = 1;
            if ~isempty(inds),
                searchables = searchables .* dict_msk;
            end
        end
        
        if isfield(priors,'flip')
            if ischar(priors.flip)
                pr.flip = deg2rad(flipmap(p));
            end
            allflip = cell2mat({parms(:).flip});
            dict_msk = zeros(size(searchables));
            inds = find( abs(allflip-pr.flip) < matchfactor*pr.flip);
            dict_msk(inds) = 1;
            if ~isempty(inds),
                
                searchables = searchables .* dict_msk;
            end
        end
        
        if isfield(priors,'cbva')
            if ischar(priors.cbva)
                pr.cbva = volmap(p);
            end
            
            allcbva = cell2mat({parms(:).cbva});
            dict_msk = zeros(size(searchables));
            inds = find( abs(allcbva-pr.cbva) < matchfactor*pr.cbva);
            dict_msk(inds) = 1;
            if ~isempty(inds),
                searchables = searchables .* dict_msk;
            end
        end
        
        if isfield(priors,'kfor')
            if ischar(priors.kfor)
                pr.kfor = kformap(p);
            end
            
            
            allkfor = cell2mat({parms(:).kfor});
            dict_msk = zeros(size(searchables));
            inds = find( abs(allkfor-pr.kfor) < matchfactor*pr.kfor);
            dict_msk(inds) = 1;
            searchables = searchables .* dict_msk;
        end
        
        %%
        
   
        
        
        inds = find(searchables);
        r_dict = dict(:,inds);
        r_parms = parms(inds);
        
        
        
        y = raw(:,p);
        % y = detrend(y);
        y = y - mean(y);
        y = y/ norm(y);
        
        
        %{
        % this is intended to normalize and mean center  each half of the
        % data separately for the cases where we're concatenating two runs
        
        y = normalize_halves(y);
       
       for n=1:size(dict,2)
           dict(:,n) = normalize_halves (dict(:,n)) ;
       end
        %}
        
        
        if ~sum(isnan(y))
            [best, R] = best_match(r_dict, y);
            bp = r_parms(best);
        else
            bp = r_parms(1);
            R = 0;
            best = 1;
        end;
        
        
        if ~isempty(best)
            
            if length(best) >= 2
                fprintf('\rWhoa ... there are %d best matches  ....   !', length(best));
            end
            
            flowmap(p) = bp(1).f;
            volmap(p) = bp(1).cbva;
            transmap(p) = bp(1).bat;
            transmap2(p) = bp(1).bat2;
            
            dispmap(p) = bp(1).Disp;
            t1map(p) = 1/bp(1).r1tis;
            Rmap(p) = R(1);
            flipmap(p) = rad2deg(bp(1).flip);
            kformap(p) = bp(1).kfor;
        end
        
        
        if Npix==1
            fprintf('\nROI analysis .... \n');
            figure
            subplot(211)
            plot(dict);
            hold on
            plot(y,'k')
            plot(y,'k*')
            title('compared to Full dictionary')
            hold off
            fprintf('\nEstimated parms:   \n');
            bp(1)
            
            subplot(212)
            plot(y,'k')
            hold on
            title('Compared to Best Entry'); pause(0.5)
            plot(r_dict(:,best))
            hold off
            drawnow
            pause(0.1)
        end
        
        r_dict=[];

    end
    
end


Rmap(isnan(Rmap)) = 0;

if Npix>1
    flowmap =  reshape(flowmap, xdim, xdim, 1);
    transmap = reshape(transmap, xdim, xdim, 1);
    transmap2 = reshape(transmap2, xdim, xdim, 1);
    dispmap = reshape(dispmap, xdim, xdim, 1);
    t1map = reshape(t1map, xdim, xdim, 1);
    Rmap = reshape(Rmap, xdim, xdim, 1);
    volmap = reshape(volmap, xdim, xdim, 1);
    flipmap = reshape(flipmap, xdim, xdim, 1);
    kformap = reshape(kformap, xdim, xdim, 1);
end

if saveflip
    save flipang.mat flipmap
end

if saveflow
    save flows.mat flowmap
end

if savebat
    save trans.mat transmap
end

if savebat2
    save trans2.mat transmap2
end

if savedisps
    save disps.mat dispmap
end

save R_flows.mat Rmap

if savekfors
    save kfors.mat kformap
end

if savet1
    save t1map.mat t1map
end

if savecbv
    save cbva.mat volmap
end

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
    lightbox((flipmap), [], 1);
    title('Flip Angle (deg)')
    
    subplot(324)
    lightbox(transmap2);
    title('BAT2(s)');
    
    subplot(325)
    lightbox(t1map);
    title('T1 map  (s)');
    
    subplot(326)
    lightbox(volmap, [], 1);
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

function junk4testing
%%
phys.kfor = 0.02;
phys.mtis0 = 1;
phys.r1tis = 0.57    ;
phys.flip = 0.91 ;
phys.f = 0.0083;
phys.cbva =      0.01 ;
phys.bat=  1.5;
phys.L = 1;
phys.Disp =      20;
phys.Ptime =     0.5;

priors.f =  0.0083;


timing_parms.PID = ones(10,1)*1.2;
timing_parms.t_tag = ones(10,1)*1.6;
timing_parms.t_adjust = ones(10,1);
timing_parms.t_aq = ones(10,1)*0.035;
timing_parms.t_delay = ones(10,1)*0.035;
timing_parms.Nlabel_group = 1;
timing_parms.order = 1;
timing_parms.isLabel = ones(10,1);
timing_parms.isLabel( [2 4 5 7 8]) = 0;


dict_phys_parms.r1tis =  1./linspace(0.5,3,5);
dict_phys_parms.flip =   deg2rad(linspace(40,90,5));
dict_phys_parms.bat = linspace(0.5, 4, 5);
dict_phys_parms.f = linspace(0,100,5) / 6000;
dict_phys_parms.cbva = 0.01; % [linspace(0.005, 0.03, 5) ];
dict_phys_parms.kfor = 0.02; % [0 0.02 0.04 0.06];
dict_phys_parms.Disp = [10 40  ];
dict_phys_parms.mtis0 = 1;

[dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);

obs = gen_signals_150521(phys, timing_parms, 1,0);

flex_search_150710(dict, parms, obs, [], priors)
%%
return
