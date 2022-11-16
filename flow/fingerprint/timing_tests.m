
% read in sequence timings
t_tag = load('t_tags.txt');
t_delay = load('t_delays.txt');
t_adjust = load('t_adjusts.txt');

% corrections - pulse sequence in the scope doesn't quite do what it's
% supposed to

t_delay = t_delay + 0.003  ;
t_adjust = t_adjust  + 0.009 ;
t_tag = t_tag - 0.002  ;

timing_parms.t_delay = t_delay;
timing_parms.t_tag = t_tag;
timing_parms.t_adjust = t_adjust;
timing_parms.Nlabel_group = 1;
timing_parms.order = 1; % use the default order (usually control-tag)

% Default order is Control-Tag
tmp = ones(size(t_adjust));
tmp(1:2:end) = 0;
timing_parms.isLabel = tmp;


if exist(fullfile(pwd,'Ntags_group.txt'))
    timing_parms.Nlabel_group = load('Ntags_group.txt');
end
if exist(fullfile(pwd,'tag_order.txt'))
    timing_parms.order = load('tag_order.txt');
end
if exist(fullfile(pwd,'labelcontrol.txt'))
    timing_parms.isLabel = load('labelcontrol.txt');
end

% Note: the pulse sequence is set up to do control in the first frame, then tag ...etc.
% if the calibration runs are negative, then we know that the order has been reversed by off-resonance
%timing_parms.order = 2; % Flip the specified order


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

obs1 = gen_signals_150521(phys_parms, timing_parms, 1, 0);

t_tag1 = timing_parms.t_tag;
t_tag2 = 0.003*floor(t_tag1 / 0.003);
timing_parms.t_tag = t_tag2;

obs2 = gen_signals_150521(phys_parms, timing_parms, 1, 0);
figure
plot(obs1(15:end)); hold on
plot(obs2(15:end),'r'); hold off


