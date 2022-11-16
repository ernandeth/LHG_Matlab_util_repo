function flex_search_160525(dict, parms, raw, msk, priors)
% function flex_search_160525(dict, parms, raw, msk, priors)
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
% flex_search_160525(dict, parms, obs, [], priors)


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
savebat = 1;
savebat2 = 1;

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
% the dictionary entries must have values within 10% of the value of the prior
matchfactor = 0.10; 

      tmp.mtis0 = [];
    tmp.f = 1;
    tmp.cbva = 1;
    tmp.bat = 1;
    tmp.bat2 = 1;
    tmp.r1tis = 1;
    tmp.flip = 1;
    tmp.Disp = 1;
    tmp.kfor = 1;

for p=1:Npix
    
    r_parms = tmp;
    
    best = 1;

    if ~mod(p,200)
        str = sprintf('echo working on pixel: %d   out of   %d   pixels', p, Npix);
        system(str);
    end
    
    if msk(p)==1 || size(raw,2)==1
        
        % searchables is a mask indicating whith entries of the dictionary
        % can be used for this pixel.  At the beginning, they are ALL
        % searchable, ie. they are all ones.
        searchables = ones(size(dict,2),1);
        
        pr = priors;
        
        % if we already the values of some parameters, then we use them 
        % to reduce the size of the dictionary.
        % if the priors are maps, then we look at each pixel.
        
        
        if isfield(priors,'f')       
            if ischar(priors.f)
                pr.f = flowmap(p);
            end
            
            % All the flow values corresponding to each entry of the
            % dictionary:
            allf = parms.f;
            
            % see which of those flow values in the dictionary are close to
            % the prior known value.  Only Those entries can be searched.
           
            inds = find( abs(allf - pr.f) < matchfactor*pr.f);
            
            % dict_msk is a "mask" indicating which dictionary entries are
            % searched for this pixel.  The mask puts zeros in the
            % "searchables" list.
            dict_msk = zeros(size(searchables));
            dict_msk(inds) = 1;
            
            if ~isempty(inds),
                searchables = searchables .* dict_msk;
            end
            
        end
        
        
        if isfield(priors,'r1tis')
            if ischar(priors.r1tis)
                pr.r1tis = 1/t1map(p);
            end
            
            allr1tis = parms.r1tis;
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
            allbat = parms.bat;
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
            allbat2 = parms.bat2;
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
            allflip = parms.flip;
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
            
            allcbva = parms.cbva;
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
            
            
            allkfor = parms.kfor;
            dict_msk = zeros(size(searchables));
            inds = find( abs(allkfor-pr.kfor) < matchfactor*pr.kfor);
            dict_msk(inds) = 1;
            searchables = searchables .* dict_msk;
        end
        
        %%
        
        
        % the indices of the entries of the dictionary that survive the reduction
        inds = find(searchables);
        
        % reduced dictionary
        r_dict = dict(:,inds); 
        reducedsize = length(inds);
        
        % reduced parameter sets corresponding to reduced dictionary 
        r_parms.mtis0 = parms.mtis0(inds);
        r_parms.f = parms.f(inds);
        r_parms.cbva = parms.cbva(inds);
        r_parms.bat = parms.bat(inds);
        r_parms.bat2 = parms.bat2(inds) ;
        r_parms.r1tis = parms.r1tis(inds);
        r_parms.flip = parms.flip(inds);
        r_parms.Disp = parms.Disp(inds);
        r_parms.kfor = parms.kfor(inds);
        
        % normalizing the raw data at each pixel:
        y = raw(:,p);
        y = y - mean(y);
        y = y/ norm(y);
        
        
        [best R]= best_match(r_dict, y);
        
        if ~isempty(best)
            
            if length(best) >= 2
                fprintf('\rWhoa ... there are %d best matches  ....   !', length(best));
                best = best(1);
            end
            
            % Extract the parameters from the reduced parameter table
            % and put them into the corresponding maps.
            flowmap(p)  = r_parms.f(best);
            volmap(p)   = r_parms.cbva(best);
            transmap(p) = r_parms.bat(best);
            transmap2(p)= r_parms.bat2(best);           
            dispmap(p)  = r_parms.Disp(best);
            t1map(p)    = 1/r_parms.r1tis(best);
            flipmap(p)  = rad2deg(r_parms.flip(best));
            kformap(p)  = r_parms.kfor(best);
            
            Rmap(p) = R;
           
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
        
        %r_dict = [];
        clear r_dict r_parms
        

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


if Npix>1
    show_parm_maps  
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
