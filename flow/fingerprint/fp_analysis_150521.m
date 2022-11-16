%using the increasing timing schedule
close all;
% matlabpool open

addpath  ~/matlab/flow/milf/

pfile = dir('P*');
%sprec1(pfile(1).name, 'fy','N');

%  read in time series and subtract
[pth nm ex] = fileparts(ls('vol*.nii'));

%spm_smooth([nm ex] , ['s_' nm ex], 8);
%nm = ['s_' nm];

[raw1 h]= read_img(nm);

aslsub(nm,1,1,h.tdim,0,0,0);
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

%matlabpool open

psdseqtime = 0.030;
nslices = 1;
doSub = 0;
doROI = 0;
search_type = 1;

fudge_mat = [];
for delay_fudge =  0; %[0:50:500]/1000
for adjust_fudge = 0; % [-15:5:15]/1000
for tag_fudge = 0 ; %[0:20:100]/1000

%
%
% read in sequence timings
t_tag = load('t_tags.txt');
t_delay = load('t_delays.txt');
t_adjust = load('t_adjusts.txt');
if exist(fullfile(pwd,'t_aqs.txt'))
	t_aq = load('t_aqs.txt');
else
    t_aq = 0.0349*ones(size(t_tag));
    t_aq = 0.0325*ones(size(t_tag));

end
%

%{ 
% using nominal schedule:
t_tag = load('nominal/t_tags.txt');
t_delay = load('nominal/t_delays.txt');
t_adjust = load('nominal/t_adjusts.txt');
if exist(fullfile(pwd,'nominal/t_aqs.txt'))
	t_aq = load('nominal/t_aqs.txt');
else
    t_aq = 0.0349*ones(size(t_tag));
end
%}

%{
t_aq =t_aq(1:60);
t_tag =t_tag(1:60);
t_delay =t_delay(1:60);
t_adjust =t_adjust(1:60);
raw = raw(1:60,:);
Nframes = 60;
%}

%{
t_aq =t_aq(2:end);
t_tag =t_tag(2:end);
t_delay =t_delay(2:end);
t_adjust =t_adjust(2:end);
raw = raw(1:end-1,:);
Nframes = Nframes-1;
%}
% corrections - pulse sequence in the scope doesn't quite do what it's
% supposed to

% t_delay = t_delay + 0.003  ;
% t_adjust = t_adjust  + 0.007 ;
% t_tag = 0.003*floor(t_tag/0.003) ;


% just trying out some other numbers...
t_delay = t_delay + delay_fudge  ;
t_adjust = t_adjust + adjust_fudge ;
t_tag = t_tag + tag_fudge ;


t_tag = t_tag(1:Nframes);
t_adjust =  t_adjust(1:Nframes);
t_delay = t_delay(1:Nframes);
t_aq = t_aq(1:Nframes);

% TR = t_adjust +  t_tag + t_delay + nslices*psdseqtime ;
% %
% duration = sum(TR)+2;

timing_parms.t_aq = t_aq;
timing_parms.t_delay = t_delay;
timing_parms.t_tag = t_tag;
timing_parms.t_adjust = t_adjust;
timing_parms.Nlabel_group = 1;
timing_parms.order = 0; % use the default order (usually control-tag)

% Default order is Control-Tag
tmp = ones(size(t_adjust));
tmp(1:2:end) = 0;
timing_parms.isLabel = tmp;


if exist(fullfile(pwd,'Ntags_group.txt'))
    timing_parms.Nlabel_group = load('Ntags_group.txt');
end
if exist(fullfile(pwd,'tag_order.txt'))
    timing_parms.order = load('tag_order.txt')
end
if exist(fullfile(pwd,'labelcontrol.txt'))
    labelcontrol =  load('labelcontrol.txt');
    timing_parms.isLabel = labelcontrol;
end
% Note: the pulse sequence is set up to do control in the first frame, then tag ...etc.
% if the calibration runs are negative, then we know that the order has been reversed by off-resonance
% timing_parms.order = 2; % Flip the specified order

%

%%  a test
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
obs = gen_signals_150521(physparms, timing_parms, 1,0);
drawnow;
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

fitmode = 'flow_t1prior2'
fitmode = 't1'
fitmode = 'flow_t1prior'
fitmode = 'flow'


switch (fitmode)
    case ('t1')
                
        dict_phys_parms.bat =  linspace(0.5,4,3);
        dict_phys_parms.bat2 = linspace(0.5,4,3) ;
        dict_phys_parms.r1tis = 1./linspace(0.5, 3.5, 20);
        dict_phys_parms.flip = deg2rad(linspace(20,60,20));
        dict_phys_parms.f = [ 0.005 0.01 ];
        dict_phys_parms.cbva = 0.01;
        dict_phys_parms.kfor = 0.02;
        dict_phys_parms.Disp = [40];
        dict_phys_parms.mtis0 = 1;
        
               
        [dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);

       
       % one dictionary for all pixels
        flex_search_150521(dict, parms, raw, msk, [], timing_parms);
        
    case ('flow')
        
        %
        dict_phys_parms.r1tis =   linspace(0.3,2,10);
        dict_phys_parms.flip =    deg2rad(linspace(30,70,10)); %  deg2rad(linspace(40,80,10));
        dict_phys_parms.bat = linspace(0.5, 4, 10);        
        dict_phys_parms.bat2 =  linspace(0.1, 4, 10) ; 
        dict_phys_parms.f =  linspace(0,150,20) / 6000; %  linspace(0,150,20) / 6000;
        dict_phys_parms.cbva =  [linspace(0.001, 0.03, 10) 0.5 0.9]; % 0.01;
        dict_phys_parms.kfor = 0.02; % [0 0.02 0.04 0.06];
        dict_phys_parms.Disp = 40;
        dict_phys_parms.mtis0 = 1;

        fprintf('\nMaking Dictionary ... ');
       
        [dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);

        fprintf('\nDictionary done');

        save -v7.3 flow_dictionary.mat dict parms
        
        tmpraw = raw(1:end,:);
        tmpdict = dict(:,1:end);
       
        %{   
        % redundant!  gen_dictionary already normalizes!
        parfor n=1:size(tmpdict,1)
            tmpdict(n,:) = normcenter(tmpdict(n,:));
        end
        %}
        
        % one dictionary for all pixels
        fprintf('\nSearching Dictionary ... ');
        flex_search_150521(tmpdict, parms, tmpraw, msk, [], timing_parms);
        
    case ('flow_t1prior')
        
        dict_phys_parms.r1tis = 'prior_t1map.mat';
        dict_phys_parms.flip = 'prior_flipang.mat';
        dict_phys_parms.bat = linspace(0.5, 4, 10);
        dict_phys_parms.bat2 =  linspace(0.5, 4, 10) ; 
        dict_phys_parms.f = linspace(0,100,20) / 6000;
        dict_phys_parms.cbva =  [linspace(0.001, 0.03, 10) 0.5 0.9]; % 0.01;
        dict_phys_parms.kfor = 0.02;
        dict_phys_parms.Disp = 40;
        dict_phys_parms.mtis0 = 1;
        
       
       % each pixel gets its own dictionary
       flex_search_150521([], [], raw, msk, dict_phys_parms, timing_parms);

    case ('flow_t1prior2')
        
        % Make a comprehensive dictionary
                
        % first search for T1, flip
        fprintf('\nStep 1: ')

        dict_phys_parms.r1tis = linspace(0.3, 3, 30);
        dict_phys_parms.flip = deg2rad(linspace(30,70,30));
        dict_phys_parms.bat =  1;
        dict_phys_parms.bat2 = 1;
        dict_phys_parms.f = [0]/6000;
        dict_phys_parms.cbva = 0.01;
        dict_phys_parms.kfor = linspace(0,1,10);
        dict_phys_parms.Disp = [40];
        dict_phys_parms.mtis0 = 1;
               
        fprintf('\nMake a small dictionary for T1 fit ... ')
        %
        [dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);
        
        fprintf('\nDictionary done');
        
        save -v7.3 t1_dictionary.mat dict parms
        load t1_dictionary.mat
        
        priors = [];
        fprintf('\nSearch for T1 fit ... ')
        flex_search_150710(dict, parms, raw, msk, priors)

        !cp t1map.mat prior_t1map.mat
        !cp flipang.mat prior_flipang.mat
        
        drawnow
        %
        
        % then search for flow and bat
        fprintf('\nStep 2: ')
        fprintf('\nMake big dictionary for Flow and everything ... ')

        dict_phys_parms.r1tis = linspace(0.3, 2, 10);
        dict_phys_parms.flip = deg2rad(linspace(30,70,10));       
        dict_phys_parms.bat =  linspace(0.2,2,10);
        dict_phys_parms.bat2 =  linspace(0.5,3,10);
        dict_phys_parms.f = linspace(0,100,20)/6000;
        dict_phys_parms.cbva = [ 2*(10.^linspace(-3,-1,10))  1];
        dict_phys_parms.kfor = 0.02;
        dict_phys_parms.Disp = [140];
        dict_phys_parms.mtis0 = 1;
        
       [dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);

        save -v7.3 whole_dictionary.mat dict parms
        load whole_dictionary.mat 

        %
        priors.r1tis ='prior_t1map.mat';
        priors.flip = 'prior_flipang.mat';
        
        %priors = [];
        %
        
        fprintf('\nSearch for Flow with T1 priors... ')
 
        % skip the first few points to downweight the importance of t1
        % info
               
        tmpraw = raw(1:end,:);
        tmpdict = dict(:,1:end);
        %{ 
        parfor n=1:size(tmpdict,1)
            tmpdict(n,:) = normcenter(tmpdict(n,:));
        end
        %}
        flex_search_150710(tmpdict, parms, tmpraw, msk, priors)

        
         fprintf('\n ...done\n ')

end

fprintf('\n del= %f tadj = %f tag = %f\n', delay_fudge, adjust_fudge, tag_fudge)
load R_flows.mat
tmp = sum(Rmap(Rmap>0))
Rfits = [Rfits ; tmp];

close
set(gcf, 'Name', sprintf('del= %f tadj = %f tag = %f', delay_fudge, adjust_fudge, tag_fudge));
drawnow

fudge_mat = [fudge_mat; delay_fudge adjust_fudge tag_fudge];

end
end
end


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
