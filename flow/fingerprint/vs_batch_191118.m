% default parameters for simulation
Nframes = 6;
dofigs = 0;
doSub= 0;

aq_parms.alpha = 0.5;
aq_parms.t_adjusts = 2*ones(Nframes,1);
aq_parms.t_tags = 0*ones(Nframes,1);
tmp = [linspace(1.3, 1.3, Nframes)']; tmp(2:2:end) = tmp(1:2:end);
aq_parms.t_delays = tmp;
aq_parms.ArtSup_delay = 0.1 ;% delay between AS pulse and acqusition
tmp =  ones(Nframes,1); %tmp(end/2+1:end) = 0;
aq_parms.doArtSup = tmp; 
tmp =  zeros(2*Nframes,1); 
tmp(3:4:end) = 1;
tmp(4:4:end) = 1;
aq_parms.labelcontrol = tmp; 

parms.f = 50/6000;
parms.Mtis0 = 1;
parms. cbva = 0.02;
parms.bat = 0.3;
parms.r1tis = 1/1.4;
parms.flip = deg2rad(90);
%%

fprintf('\nsingle VSS experiment:\n')
S = gen_signals_vs_191105(parms, aq_parms, 1,0);
deltaS = mean(S(4:2:end) - S(3:2:end))
S0 = mean(S(4:2:end))
pct_signal = 100*deltaS/S0

%%
fprintf('\nsingle VSI experiment:\n')
aq_parms.alpha = 0.9;
S = gen_signals_vs_191105(parms, aq_parms, 1,0);
deltaS = mean(S(4:2:end) - S(3:2:end))
S0 = mean(S(4:2:end))
pct_signal = 100*deltaS/S0
%%
fprintf('\ndual VSS experiment:\n')
aq_parms.alpha = 0.5;
aq_parms.t_tags = 1.1*ones(Nframes,1);
aq_parms.t_delays(:) = 0.8;

S = gen_signals_vs_191105(parms, aq_parms, 1,0);
deltaS = mean(S(4:2:end) - S(3:2:end))
S0 = mean(S(4:2:end))
pct_signal = 100*deltaS/S0
%%
fprintf('\ndual VSI experiment:\n')
aq_parms.alpha = 0.9;
aq_parms.t_tags = 1.8*ones(Nframes,1);
aq_parms.t_delays(:) = 0.5;

S = gen_signals_vs_191105(parms, aq_parms, 1,0);
deltaS = mean(S(4:2:end) - S(3:2:end))
S0 = mean(S(4:2:end))
pct_signal = 100*deltaS/S0
%%

fprintf('\nOptimal t_delay and t_tag for dual saturations.\n')
Nframes = 80;
dofigs = 0;
doSub= 0;

aq_parms.alpha = 0.5;
aq_parms.t_adjusts = 2*ones(Nframes,1);
tmp = [linspace(0.5, 3, Nframes)']; tmp(2:2:end) = tmp(1:2:end);
aq_parms.t_delays = tmp;
aq_parms.ArtSup_delay = 0.1 ;% delay between AS pulse and acqusition
tmp =  ones(Nframes,1); %tmp(end/2+1:end) = 0;
aq_parms.doArtSup = tmp; 
tmp =  zeros(2*Nframes,1); 
tmp(3:4:end) = 1;
tmp(4:4:end) = 1;
aq_parms.labelcontrol = tmp; 

all_pct = [];
t_tags = [0.2:0.1:2.5]';
for c=1:length(t_tags)
    aq_parms.t_tags = t_tags(c)*ones(Nframes,1);
    S = gen_signals_vs_191105(parms, aq_parms, 0,0);
    deltaS = (S(4:2:end) - S(3:2:end));  % skip the first pair
    S0 =   (S(4:2:end));    % skip the first pair
    pct_signal = abs(100*deltaS./S0);
    all_pct = [all_pct; pct_signal ];
end
imagesc(all_pct)
caxis([-1 1]*0.8) 
colorbar
xlabel('post label delay')
ylabel('label duration')
msig = max(all_pct(:))
c=find(all_pct==msig);
[maxdur maxdel] = ind2sub(size(all_pct),c)

fprintf('\nOptimal t_delay = %f and t_tag =%f for dual saturations.\n', ...
    aq_parms.t_delays(maxdel+2), t_tags(maxdur))
fprintf('\n percent signal = %f \n', msig)

%%

fprintf('\nOptimal t_delay for dual inversions.\n')
aq_parms.alpha = 0.9;

Nframes = 80;
dofigs = 0;
doSub= 0;

aq_parms.t_adjusts = 2*ones(Nframes,1);
tmp = [linspace(0.6, 3, Nframes)']; tmp(2:2:end) = tmp(1:2:end);
aq_parms.t_delays = tmp;
aq_parms.ArtSup_delay = 0.1 ;% delay between AS pulse and acqusition
tmp =  ones(Nframes,1); %tmp(end/2+1:end) = 0;
aq_parms.doArtSup = tmp; 
tmp =  zeros(2*Nframes,1); 
tmp(3:4:end) = 1;
tmp(4:4:end) = 1;
aq_parms.labelcontrol = tmp; 

all_pct = [];
t_tags = [0.4:0.1:2]';
for c=1:length(t_tags)
    aq_parms.t_tags = t_tags(c)*ones(Nframes,1);
    S = gen_signals_vs_191105(parms, aq_parms, 0,0);
    deltaS = (S(4:2:end) - S(3:2:end));   % skip first pair
    S0 = (S(4:2:end));  % skip first pair
    pct_signal = (100*deltaS./S0);
    all_pct = [all_pct; pct_signal ];
end
imagesc(all_pct)
caxis([-1 1]*0.8) 
colorbar
xlabel('post label delay')
ylabel('label duration')
msig = max(abs(all_pct(:)))
c=find(abs(all_pct)==msig);
[maxdur maxdel] = ind2sub(size(all_pct),c)

fprintf('\nOptimal t_delay = %f and t_tag =%f for dual inversions.\n', ...
    aq_parms.t_delays(maxdel+2), t_tags(maxdur))
fprintf('\n percent signal = %f \n', all_pct(c))


%% Check sensitivity to perfusion
parms.f = 50/6000;
obs0 = gen_signals_vs_191105(parms, aq_parms, 0,0)
parms.f = 0;
obs1 = gen_signals_vs_191105(parms, aq_parms, 0,0)

plot(obs1-obs0)


%% make some fingerprints
aq_parms.t_tags = 1.2*ones(Nframes,1);
tmp = 0.2 + rand(Nframes,1)*2;
aq_parms.t_delays = tmp;

parms.f = 50/6000;
obs0 = gen_signals_vs_191105(parms, aq_parms, 0,0)
parms.f = 0;
obs1 = gen_signals_vs_191105(parms, aq_parms, 0,0)

plot(obs1-obs0)


