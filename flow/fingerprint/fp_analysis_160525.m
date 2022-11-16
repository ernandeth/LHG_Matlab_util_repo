%using the increasing timing schedule
close all;
% matlabpool open

% startup_hernan

pfile = dir('P*');
%sprec1(pfile(1).name, 'fy','N');

%  read in time series and subtract
[pth nm ex] = fileparts(ls('vol*.nii'));  % <---
vfiles = dir('vol*.nii');
nm = vfiles(1).name;

%spm_smooth([nm ex] , ['s_' nm ex], 8);
%nm = ['s_' nm];

[raw1 h]= read_img(nm);

aslsub(nm(1:end-4),1,1,h.tdim,0,0,0);
raw1 = read_img('sub.img');


[tmp hdr] = read_img('mean_con.img');
msk = tmp;
msk(msk<3000) = 0;
msk(msk>0) = 1;
lightbox(reshape(msk, hdr.xdim, hdr.ydim, hdr.zdim));
drawnow

raw = read_img(nm);
Nframes = size(raw,1);
Rfits = [];
Npix = size(raw, 2);

% this bit will regress out the average signal outside the mask assuming
% it's the fat/off-resonance noise
%{
noisevec = 0;
for p=1:Npix
    if msk(p)==0 & ~isnan(msk(p))
        noisevec = noisevec + raw(:,p);
    end
end

noisevec = noisevec- mean(noisevec);
raw = raw - noisevec*pinv(noisevec)*raw;
        
%}

%matlabpool open

psdseqtime = 0.030;
nslices = 1;
doSub = 0;
doROI = 0;
search_type = 1;

fudge_mat = [];
% for delay_fudge =  0; %[0:50:500]/1000
% for adjust_fudge = 0; % [-15:5:15]/1000
% for tag_fudge = 0 ; %[0:20:100]/1000

%
%
% read in sequence timings
% t_tag = load('t_tags.txt');
% t_delay = load('t_delays.txt');
% t_adjust = load('t_adjusts.txt');
% 
% if exist(fullfile(pwd,'t_aqs.txt'))
%     t_aq = load('t_aqs.txt');
% else
%     t_aq = 0.0329*ones(size(t_tag));
% end
%

%
% using nominal schedule:
t_tag = load('nominal/t_tags.txt');
t_delay = load('nominal/t_delays.txt');
t_adjust = load('nominal/t_adjusts.txt');
if exist(fullfile(pwd,'nominal/t_aqs.txt'))
	t_aq = load('nominal/t_aqs.txt');
else
    t_aq = 0.0329*ones(size(t_tag));
end
%

% corrections - pulse sequence in the scope doesn't quite do what it's
% supposed to :
t_tag = 0.003*floor(t_tag/0.003) ;  % just roundoff
t_delay = t_delay + 0.005  ; % to the center of the RF pulse
t_aq = t_aq - 0.005; % from center of flip to end of crusher
t_adjust = t_adjust  + 0.003 ; % from crusher to tag

% just trying out some other numbers...
% t_delay = t_delay + delay_fudge  ;
% t_adjust = t_adjust + adjust_fudge ;
% t_tag = t_tag + tag_fudge ;


%% try to reduce data size 
Nframes = size(raw, 1);
% Nframes = Nframes/2;              % <--- * half data analysis
% Nframes = Nframes/4;              % <--- *  quarter data analysis

raw = raw(1:Nframes, :);       

t_tag = t_tag(1:Nframes);
t_adjust =  t_adjust(1:Nframes);
t_delay = t_delay(1:Nframes);
t_aq = t_aq(1:Nframes);
%%


% TR = t_adjust +  t_tag + t_delay + nslices*psdseqtime ;
% %
% duration = sum(TR)+2;

timing_parms.t_aq = t_aq(:);
timing_parms.t_delay = t_delay(:);
timing_parms.t_tag = t_tag(:);
timing_parms.t_adjust = t_adjust(:);
timing_parms.Nlabel_group = 1;
timing_parms.order = 0; % use the default order (usually control-tag)

% Default order is Control-Tag
tmp = ones(size(t_adjust));
tmp(1:2:end) = 0;
timing_parms.isLabel = tmp(:);


if exist(fullfile(pwd,'nominal/Ntags_group.txt'))
    timing_parms.Nlabel_group = load('nominal/Ntags_group.txt');
end
if exist(fullfile(pwd,'tag_order.txt'))
    timing_parms.order = load('tag_order.txt');
end
if exist(fullfile(pwd,'nominal/labelcontrol.txt'))
    labelcontrol =  load('nominal/labelcontrol.txt');
    timing_parms.isLabel = labelcontrol(:);
end
% Note: the pulse sequence is set up to do control in the first frame, then tag ...etc.
% if the calibration runs are negative, then we know that the order has been reversed by off-resonance
% timing_parms.order = 2; % Flip the specified order

%

%  a test
%
physparms.f=         60 /6000;
physparms.mtis0 =     1 ;
physparms.cbva =      0.02 ;
physparms.bat =   1.2 ;
physparms.bat2 =   0.5 ;

physparms.kfor =      1e-2 ;
physparms.r1tis =     1/1.4  ;
physparms.flip =      45*pi/180 ; % flip angle in radians
physparms.Disp =      20;

%obs = gen_signals_150521(physparms, timing_parms, 1,0);
obs = gen_signals_160426(physparms, timing_parms, 1,0);

drawnow;

% return
%
%%

%{
sprintf('\nDoing an ROI ....');
raw = load ('GM_rawtdata.dat');
%}

Npix = size(raw,2);
raw = single(raw);

n=0;
converged=0;
oldfit = 0;

if batchmode==0
    fitmode = 'flow'
    fitmode = 't1'
    fitmode = 'flow_t1prior'
    fitmode = 'flow_t1prior2'
    nominalFlip = 60;
end



switch (fitmode)
    case ('t1')
        
        dict_phys_parms.bat =  linspace(0.5,4,3);
        dict_phys_parms.bat2 = linspace(0.5,4,3) ;
        dict_phys_parms.r1tis = 1./linspace(0.5, 3.5, 20);
        dict_phys_parms.flip = deg2rad(linspace(20,60,20));
        dict_phys_parms.f = [ 0.005 0.01 ];
        dict_phys_parms.cbva = 0.01;
        dict_phys_parms.kfor = 0.0;
        dict_phys_parms.Disp = [40];
        dict_phys_parms.mtis0 = 1;
        
        % [dict, parms] = gen_flex_dictionary_150526 (timing_parms, dict_phys_parms);
        [dict, parms] = gen_flex_dictionary_160525 (timing_parms, dict_phys_parms);
        
        
        % one dictionary for all pixels
        % flex_search_150521(dict, parms, raw, msk, [], timing_parms);
        flex_search_160525(dict, parms, raw, msk, []);
        
    case ('flow')
        
        
        dict_phys_parms.r1tis = linspace(0.35, 1.5, 10);
        dict_phys_parms.flip = deg2rad(linspace(0.6*nominalFlip, 1.4*nominalFlip,5));
        dict_phys_parms.bat =  linspace(0.5,3,10);
        dict_phys_parms.bat2 =  linspace(0.5,3,10);
        dict_phys_parms.f = linspace(0,120,10)/6000;
        dict_phys_parms.cbva = [ linspace(0,0.02,5)  1];
        dict_phys_parms.kfor =  0.02;
        dict_phys_parms.Disp = [40];
        dict_phys_parms.mtis0 = 1;
        
        fprintf('\nMaking dictionary ... ');
        
        [dict, parms] = gen_flex_dictionary_160525 (timing_parms, dict_phys_parms);
        
        save -v7.3 flow_dictionary.mat dict parms
        
        load flow_dictionary.mat
        
        tmpraw = raw(1:end,:);
        tmpdict = dict(:,1:end);
        
        
        
        % one dictionary for all pixels
        fprintf('\nSearching Dictionary ... ');
        % flex_search_150521(dict, parms, raw, msk, [], timing_parms);
        flex_search_160525(dict, parms, raw, msk, []);
        
    case ('flow_t1prior')        
        !mkdir results

        % Make a smoothing kernel:
        F = 0.5*ones(3);
        F(2,2) = 1;
        F = F/sum(F(:));
        
        fprintf('\nSearching for T1, B1 ...');
        
        dict_phys_parms.bat =  linspace(0.5,3,3);
        dict_phys_parms.bat2 = linspace(0.5,3,3) ;
        dict_phys_parms.r1tis = 1./linspace(0.5, 3.5, 40);
        dict_phys_parms.flip = deg2rad(linspace(0.5*nominalFlip, 1.5*nominalFlip,40));
        dict_phys_parms.f = [ 0.005 0.01 ];
        dict_phys_parms.cbva = 0.01;
        dict_phys_parms.kfor = linspace(0.05, 0.2, 20);
        dict_phys_parms.Disp = [40];
        dict_phys_parms.mtis0 = 1;
        
        calcT1B1 = 1;
        
        if calcT1B1
            % if we are not using precomputed T1 and B1 from alternative acquistion
            
%             [dict, parms] = gen_flex_dictionary_160525 (timing_parms, dict_phys_parms);
%             flex_search_160525(dict, parms, raw, msk, []);

                    
            % flex_search_150521([], [], raw, msk, dict_phys_parms, timing_parms);
            % one dictionary for all pixels
            [dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);
            flex_search_150521(dict, parms, raw, msk, [], timing_parms);
            
        end
        
        fprintf('\nSearching for BAT, CBV ...');
                
        % First save and smooth the prior maps
        !cp t1map.mat results/
        !cp flipang.mat results/        
        !cp kfors.mat results/        

        
        %load('t1map.mat');  t1map = filter2(F,t1map); save t1map.mat t1map
        load('flipang.mat');  flipmap = filter2(F,flipmap); save sm_flipang.mat flipmap
        load('kfors.mat');  kformap = filter2(F,kformap); save sm_kfors.mat kformap

        dict_phys_parms.r1tis = 't1map.mat';
        dict_phys_parms.flip = 'sm_flipang.mat';
        dict_phys_parms.kfor = 'sm_kfors.mat';
        % dict_phys_parms.bat2 =  linspace(0.5, 2, 20) ;
        dict_phys_parms.f = linspace(10,100,3) / 6000;
        dict_phys_parms.bat = linspace(0.5, 3, 30);
        dict_phys_parms.cbva =  [linspace(0.001, 0.03, 20) 0.5 0.9]; % 0.01;

        % each pixel gets its own dictionary
        flex_search_150521([], [], raw, msk, dict_phys_parms, timing_parms);
        
        
        fprintf('\nSearching for BAT2, CBF ...');

        % First save the prior maps tp a safe place
        !cp trans.mat results/
        !cp cbva.mat results/
        % !cp trans2.mat results/
        
        % now smooth them for the next step
        load('trans.mat');  transmap = filter2(F,transmap); save sm_trans.mat transmap
        load('cbva.mat');  volmap = filter2(F,volmap); save sm_cbva.mat volmap
        %load('trans2.mat');  transmap2 = filter2(F,transmap2); save trans2.mat transmap2
   
        dict_phys_parms.cbva = 'sm_cbva.mat'; 
        dict_phys_parms.bat = 'sm_trans.mat';        
        %dict_phys_parms.bat2 = 'trans2.mat';        
        %dict_phys_parms.cbva =  [linspace(0.001, 0.03, 20) 0.5 0.9]; % 0.01;
        dict_phys_parms.bat2 =  linspace(0.5, 3, 30) ;
        dict_phys_parms.f = linspace(0,150,30) / 6000;
        
        flex_search_150521([], [], raw, msk, dict_phys_parms, timing_parms);
        
        % now save the new results to a safe place
        !cp flows.mat results/
        %!cp cbva.mat results/
        !cp trans2.mat results/
        !cp R_flows.mat results/
        !cp kfors.mat results/

        cd results
        show_parm_maps
        cd ..
        %{
        
        priors.r1tis = 'prior_t1map.mat';
        priors.flip = 'prior_flipang.mat';
       
        fprintf('\nMake big dictionary for all parms ... ')
        tic
        [dict, parms] = gen_flex_dictionary_160525 (timing_parms, dict_phys_parms);
        toc
        
        save -v7.3  dictionary_t1prior.mat dict parms
        
        fprintf('\nsearching dictionary with priors ...')
        tic
        flex_search_160525(dict, parms,raw, msk, priors);
        toc
        %}
        
    case ('flow_t1prior2')
        
       % This case is a second pass using the previously computed maps as priors. 
        !mkdir results_pass2
        
        % use first pass prior maps
        !cp results/t1map.mat results_pass2/
        !cp results/flipang.mat results_pass2/
        !cp results/kfors.mat results_pass2/

        !cp results/t1map.mat .
        !cp results/flipang.mat .
        !cp results/kfors.mat .
        
        % default search space
        dict_phys_parms.bat =  linspace(0.5,3,3);
        dict_phys_parms.bat2 = linspace(0.5,3,3) ;
        dict_phys_parms.r1tis = 1./linspace(0.5, 3.5, 30);
        dict_phys_parms.flip = deg2rad(linspace(20,60,30));
        dict_phys_parms.f = [ 0.005 0.01 ];
        dict_phys_parms.cbva = 0.01;
        dict_phys_parms.kfor = 0.1;
        dict_phys_parms.Disp = [40];
        dict_phys_parms.mtis0 = 1;
        
        
        fprintf('\nSearching for BAT, CBV - second pass...');
        fprintf('\n(using the BAT2, Flow, T1 and B1 from the first pass)...');
        
        % First save and smooth the prior maps
        % Make a smoothing kernel:
        F = 0.5*ones(3);
        F(2,2) = 1;
        F = F/sum(F(:));
                
        load('results/flipang.mat');  flipmap = filter2(F,flipmap); save sm_flipang.mat flipmap
        load('results/flows.mat');  flowmap = filter2(F,flowmap); save sm_flows.mat flowmap
        load('results/trans2.mat');  transmap = filter2(F,transmap2); save sm_trans2.mat transmap2
        load('kfors.mat');  kformap = filter2(F,kformap); save sm_kfors.mat kformap
        
        dict_phys_parms.r1tis = 't1map.mat';
        dict_phys_parms.flip =  'sm_flipang.mat';
        dict_phys_parms.bat2 =  'sm_trans2.mat' ;
        dict_phys_parms.f =     'sm_flows.mat';
        dict_phys_parms.kfor = 'sm_kfors.mat';     
        dict_phys_parms.bat = linspace(0.5, 3, 40);
        dict_phys_parms.cbva =  [linspace(0.001, 0.03, 40) 0.5 0.9]; % 0.01;

        % each pixel gets its own dictionary
        flex_search_150521([], [], raw, msk, dict_phys_parms, timing_parms);
        
        show_parm_maps
        
        fprintf('\nSearching for BAT2, CBF ...');

        % First save and smooth the prior maps
        !cp trans.mat results_pass2/
        !cp cbva.mat results_pass2/
               
        load('trans.mat');  transmap = filter2(F,transmap); save sm_trans.mat transmap
        load('cbva.mat');  volmap = filter2(F,volmap); save sm_cbva.mat volmap
   
        dict_phys_parms.cbva = 'sm_cbva.mat'; 
        dict_phys_parms.bat = 'sm_trans.mat';        
        dict_phys_parms.bat2 =  linspace(0.5, 3, 40) ;
        dict_phys_parms.f = linspace(0,150,40) / 6000;
        
        flex_search_150521([], [], raw, msk, dict_phys_parms, timing_parms);
        
        !cp flows.mat results_pass2/
        %!cp cbva.mat results/
        !cp trans2.mat results_pass2/
        !cp R_flows.mat results_pass2/
        
        cd results_pass2
        show_parm_maps
        cd ..
                
        fprintf('\n ...done\n ')

        
end

%fprintf('\n del= %f tadj = %f tag = %f\n', delay_fudge, adjust_fudge, tag_fudge)
% load R_flows.mat
% tmp = sum(Rmap(Rmap>0))
% Rfits = [Rfits ; tmp];
% 
% close
%set(gcf, 'Name', sprintf('del= %f tadj = %f tag = %f', delay_fudge, adjust_fudge, tag_fudge));
drawnow

%fudge_mat = [fudge_mat; delay_fudge adjust_fudge tag_fudge];
%
% end
% end
% end


%plot(Rfits)

%{
% Just a quick test to make sure everything works as expected
% set up default values for fit
phys_parms.f =         0.008;
phys_parms.mtis0 =     1 ;
phys_parms.cbva =       0.01;
phys_parms.bat=    1 ;
phys_parms.kfor =     0.02 ;
phys_parms.r1tis =     1/1.4  ;
phys_parms.flip =      50*pi/180 ; % flip angle in radians
phys_parms.L = 1;
phys_parms.Disp =      80;  % used to be 30
phys_parms.Ptime =     0.5;

obs = gen_signals_150521(phys_parms, timing_parms, 1, 0);
drawnow
if doSub
    obs = differencer(obs,4);
end
obs = obs - mean(obs);
obs = obs / norm(obs);
%}
