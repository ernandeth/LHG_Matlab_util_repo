% using nominal schedule:
t_tag = load('nominal/t_tags.txt');
t_delay = load('nominal/t_delays.txt');
t_adjust = load('nominal/t_adjusts.txt');
if exist(fullfile(pwd,'nominal/t_aqs.txt'))
	t_aq = load('nominal/t_aqs.txt');
else
    t_aq = 0.0329*ones(size(t_tag));
end
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
%



% corrections - pulse sequence in the scope doesn't quite do what it's
% supposed to :
t_tag = 0.003*floor(t_tag/0.003) ;  % just roundoff
t_delay = t_delay + 0.005  ; % to the center of the RF pulse
t_aq = t_aq - 0.005; % from center of flip to end of crusher
t_adjust = t_adjust  + 0.003 ; % from crusher to tag


timing_parms.t_aq = t_aq(:);
timing_parms.t_delay = t_delay(:);
timing_parms.t_tag = t_tag(:);
timing_parms.t_adjust = t_adjust(:);
timing_parms.Nlabel_group = 1;
timing_parms.order = 0; % use the default order (usually control-tag)


% Note: the pulse sequence is set up to do control in the first frame, then tag ...etc.
% if the calibration runs are negative, then we know that the order has been reversed by off-resonance
% timing_parms.order = 2; % Flip the specified order

% default physio parms


physparms.f=         60 /6000;
physparms.mtis0 =     1 ;
physparms.cbva =      0.02 ;
physparms.bat =   1.2 ;
physparms.bat2 =  1.2 ;
physparms.kfor =      1e-2 ;
physparms.r1tis =     1/1.4  ;
physparms.flip =      60*pi/180 ; % flip angle in radians
physparms.Disp =      20;

default_physparms = physparms;

%%
f = [0:5:100]/6000;
dSig = zeros(size(f));
physparms = default_physparms;
for n=1:length(f);
    physparms.f = f(n);
    tmp_phys1 = physparms;
    tmp_phys2 = physparms;
    tmp_phys1.f = physparms.f + 0.025*f(n);
    tmp_phys2.f = physparms.f - 0.025*f(n);
    
    s0 = gen_signals_160426(physparms , timing_parms, 0, 0);
    s1 = gen_signals_160426(tmp_phys1 , timing_parms, 0, 0);
    s2 = gen_signals_160426(tmp_phys2 , timing_parms, 0, 0);
    
    dSig(n) = mean(abs(s2-s1)./s0)/(0.05*f(n));
    dSig(n) = norm( (s1-s2)) / (0.05*f(n)) / norm(s0);
end
figure(432), subplot(4,2,1), 
plot(f*6000,dSig)
title('dS/df')

%%
cbva = [0:5:100]*1e-3
dSig = zeros(size(cbva));
physparms = default_physparms;
for n=1:length(cbva);
    physparms.cbva = cbva(n);
    tmp_phys1 = physparms;
    tmp_phys2 = physparms;
    tmp_phys1.cbva = physparms.cbva + 0.025*cbva(n);
    tmp_phys2.cbva = physparms.cbva - 0.025*cbva(n);
    
    s0 = gen_signals_160426(physparms , timing_parms, 0, 0);
    s1 = gen_signals_160426(tmp_phys1 , timing_parms, 0, 0);
    s2 = gen_signals_160426(tmp_phys2 , timing_parms, 0, 0);
    
    dSig(n) = mean(abs(s2-s1)./s0)/(0.05*cbva(n));
    
    dSig(n) = norm( (s1-s2)/ (0.05*cbva(n))) /norm(s0);

end
figure(432), subplot(4,2,2), 
plot(cbva,dSig)
title('dS/dCBV')

%%
r1tis = linspace(0.1, 2, 21);
dSig = zeros(size(r1tis)) ;
physparms = default_physparms;   
for n=1:length(r1tis);
    physparms.r1tis = r1tis(n);
    tmp_phys1 = physparms;
    tmp_phys2 = physparms;
    tmp_phys1.r1tis = physparms.r1tis + 0.025*r1tis(n);
    tmp_phys2.r1tis = physparms.r1tis - 0.025*r1tis(n);
    
    s0 = gen_signals_160426(physparms , timing_parms, 0, 0);
    s1 = gen_signals_160426(tmp_phys1 , timing_parms, 0, 0);
    s2 = gen_signals_160426(tmp_phys2 , timing_parms, 0, 0);
    
    dSig(n) = mean(abs(s2-s1)./s0)/(0.05*r1tis(n));

    dSig(n) = norm( (s1-s2)) / (0.05*r1tis(n)) / norm(s0);
end
figure(432), subplot(4,2,3)
plot(r1tis,dSig)
title('dS/dR1')

%%
flip = linspace(0.1, pi/2, 21);
dSig = zeros(size(flip));
physparms = default_physparms;    
for n=1:length(flip);
    physparms.flip = flip(n);
    tmp_phys1 = physparms;
    tmp_phys2 = physparms;
    tmp_phys1.flip = physparms.flip + 0.025*flip(n);
    tmp_phys2.flip = physparms.flip - 0.025*flip(n);
    
    s0 = gen_signals_160426(physparms , timing_parms, 0, 0);
    s1 = gen_signals_160426(tmp_phys1 , timing_parms, 0, 0);
    s2 = gen_signals_160426(tmp_phys2 , timing_parms, 0, 0);
    
    dSig(n) = mean(abs(s2-s1)./s0)/(0.05*flip(n));

    dSig(n) = norm( (s1-s2)) / (0.05*flip(n)) / norm(s0);

end
figure(432), subplot(4,2,4)
plot(flip*180/pi,dSig)
title('dS/dflip')
%%
bat1 = linspace(0.1, 3, 21);
dSig = zeros(size(bat1));
physparms = default_physparms;    
for n=1:length(bat1);
    physparms.bat = bat1(n);
    tmp_phys1 = physparms;
    tmp_phys2 = physparms;
    tmp_phys1.bat = physparms.bat + 0.025 * bat1(n);
    tmp_phys2.bat = physparms.bat - 0.025 * bat1(n);
    
    s0 = gen_signals_160426(physparms , timing_parms, 0, 0);
    s1 = gen_signals_160426(tmp_phys1 , timing_parms, 0, 0);
    s2 = gen_signals_160426(tmp_phys2 , timing_parms, 0, 0);
    
    dSig(n) = mean(abs(s2-s1)./s0)/(0.05*bat1(n));    
    
    dSig(n) = norm( (s1-s2)) / (0.05*bat1(n)) /norm(s0);

end
figure(432), subplot(4,2,5)
plot(bat1,dSig)
title('dS/dBAT1')
 
%%
bat2 = linspace(0.1, 3, 21);
dSig = zeros(size(bat2));
physparms = default_physparms;    
for n=1:length(bat2);
    physparms.bat2 = bat2(n);
    tmp_phys1 = physparms;
    tmp_phys2 = physparms;
    tmp_phys1.bat2 = physparms.bat2 + 0.025*bat2(n);
    tmp_phys2.bat2 = physparms.bat2 - 0.025*bat2(n);
    
    s0 = gen_signals_160426(physparms , timing_parms, 0, 0);
    s1 = gen_signals_160426(tmp_phys1 , timing_parms, 0, 0);
    s2 = gen_signals_160426(tmp_phys2 , timing_parms, 0, 0);
    
    dSig(n) = mean(abs(s2-s1)./s0)/(0.05*bat2(n));

    dSig(n) = norm( (s1-s2)) / (0.05*bat2(n)) / norm(s0);
end
figure(432), subplot(4,2,6)
plot(bat2,dSig)
title('dS/dBAT2')

%%
kfor = linspace(0, 2, 21);
dSig = zeros(size(kfor));
physparms = default_physparms;    
for n=1:length(kfor);
    physparms.kfor = kfor(n);
    tmp_phys1 = physparms;
    tmp_phys2 = physparms;
    tmp_phys1.kfor = physparms.kfor + 0.025*kfor(n);
    tmp_phys2.kfor = physparms.kfor - 0.025*kfor(n);
    
    s0 = gen_signals_160426(physparms , timing_parms, 0, 0);
    s1 = gen_signals_160426(tmp_phys1 , timing_parms, 0, 0);
    s2 = gen_signals_160426(tmp_phys2 , timing_parms, 0, 0);
    
    dSig(n) = mean(abs(s2-s1)./s0)/(0.05*kfor(n));

    dSig(n) = norm( (s1-s2) ) / (0.05*kfor(n)) /norm(s0);
end
figure(432), subplot(4,2,7)
plot(kfor,dSig)
title('dS/dMT')


%%
Disp = linspace(0, 100, 21);
dSig = zeros(size(Disp));
physparms = default_physparms;    
for n=1:length(Disp);
    physparms.Disp = Disp(n);
    tmp_phys1 = physparms;
    tmp_phys2 = physparms;
    tmp_phys1.Disp = physparms.Disp + 0.025*Disp(n);
    tmp_phys2.Disp = physparms.Disp - 0.025*Disp(n);
    
    s0 = gen_signals_160426(physparms , timing_parms, 0, 0);
    s1 = gen_signals_160426(tmp_phys1 , timing_parms, 0, 0);
    s2 = gen_signals_160426(tmp_phys2 , timing_parms, 0, 0);
    
    dSig(n) = mean(abs(s2-s1)./s0)/(0.05*Disp(n));

    dSig(n) = norm( (s1-s2)) / (0.05*Disp(n)) /norm(s0);
end
figure(432), subplot(4,2,8)
plot(Disp,dSig)
title('dS/dDisp')
