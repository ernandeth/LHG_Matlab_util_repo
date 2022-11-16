function [G, xcmat] =dictionary_correlation_150731 (timing_parms)
% function [G, xcmat] = dictionary_correlation_150731(timing_parms)
%
% this function  takes in a timing paramaters STRUCTURE
% for the ASL fingerprint experiment
% and uses it to create a dictionary of ASL signals.
% it returns the correlation among the dictionary entries
%
%
% it returns the norm of the correlation matrix of the dictionary
%
close all

if nargin==0
    timing_parms = read_timing_parms
end

%
% Combinations of  Physiological parameters for parameters:
dict_phys_parms.r1tis = 1/1.3; % 1./linspace(0.5, 3, 3);
dict_phys_parms.flip =   deg2rad(60); % deg2rad([50, 60, 70]);
dict_phys_parms.bat =  linspace(0.5, 4, 5);
dict_phys_parms.bat2 =  linspace(0.5, 4, 5);
dict_phys_parms.f =  linspace(0,100,5) / 6000;
dict_phys_parms.cbva = [linspace(0, 0.05, 5)];% , linspace(0.035,1,5)];
dict_phys_parms.kfor = 0.02;
dict_phys_parms.Disp = 20;
dict_phys_parms.mtis0 = 1;

% For debugging:
if 0
    phys.r1tis = 1;
    phys.flip =  deg2rad(60);
    phys.bat =  1;
    phys.f =  0.01;
    phys.cbva = 0.01;
    phys.kfor = 0.02;
    phys.Disp = 20;
    phys.mtis0 = 1;
    entry = gen_signals_160426( phys , timing_parms, 1, 0);
    %plot(timing_parms.t_tag); drawnow
end

[dict, parms] = gen_flex_dictionary_160525 (timing_parms, dict_phys_parms);

sigma2 = sum(var(dict'));

xc = abs(corrcoef(dict'));
xcmat = xc;

figure
for n=1:length(xcmat), xcmat(n,n) = nan; end
subplot(121); imagesc(xcmat); axis square; colormap gray; colorbar
title('Correlations Among Entries')

subplot(122); hist(xcmat(:), 100); axis square
title('Distribution of correlations')
xlabel('Correlation Coefficient')
ylabel('Number of ocurrences')

% I only care about half of the matrix
for n=1:size(xc,1)
    xc(n, 1:n) = 0;
end

% penalize entries greater than 0.9 correlation
xc(xc>0.9) = 1.5*xc(xc>0.9);

xc = (xc(:));
% G = mean(xc);
% I'm adding a penalty for having too little variance.
G = mean(xc) + 0.5/sigma2;
return
    

function timing_parms = read_timing_parms
% Loading timing parameters from files.

isLabel = load('labelcontrol.txt');
t_tag = load('t_tags.txt');
t_delay = load('t_delays.txt');
t_adjust = load('t_adjusts.txt');
if exist(fullfile(pwd,'t_aqs.txt'))
	t_aq = load('nominal/t_aqs.txt');
else
    t_aq = 0.0329*ones(size(t_tag));
end

% corrections - pulse sequence in the scope doesn't quite do what it's
% supposed to :
t_tag = 0.003*floor(t_tag/0.003) ;  % just roundoff
t_delay = t_delay + 0.005  ; % to the center of the RF pulse
t_aq = t_aq - 0.005; % from center of flip to end of crusher
t_adjust = t_adjust  + 0.003 ; % from crusher to tag
% 

% reduce data size to Nframes:
Nframes = 500;

t_tag = t_tag(1:Nframes);
t_adjust =  t_adjust(1:Nframes);
t_delay = t_delay(1:Nframes);
t_aq = t_aq(1:Nframes);
isLabel = isLabel(1:500);

timing_parms.isLabel = isLabel(:);
timing_parms.t_aq = t_aq(:);
timing_parms.t_delay = t_delay(:);
timing_parms.t_tag = t_tag(:);
timing_parms.t_adjust = t_adjust(:);
timing_parms.Nlabel_group = 1;
timing_parms.order = 0; % use the default order (usually control-tag)
