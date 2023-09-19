Nframes = 4;
TR = 4;
doSub = 1;
dofigs = 1;
aqwindow = 0.6;

parms.f=         60 /6000;
parms.mtis0 =     1 ;
parms.cbva =       0.01 ;
parms.bat =  0.1; % 0.8 ;
parms.bat2 = 1.2;
parms.kfor =      0 ;
parms.r1tis =     1/1.4  ;

parms.flip =      40*pi/180 ; % flip angle in radians
parms.Disp =      40;

 
timing_parms.t_aq = aqwindow;
timing_parms.t_delay  = 0.6 * ones(Nframes,1);
timing_parms.t_tag = 3.6*ones(Nframes,1);
timing_parms.t_adjust = TR - timing_parms.t_delay - timing_parms.t_aq;
timing_parms.order = 1;  % control-tag
timing_parms.Nlabel_group = 1;  % Number of tags in a row  (added this on 5/7/15)

tmp = ones(Nframes,1);
tmp(1:2:end) = 0;
timing_parms.isLabel = tmp';
        
obs = gen_signals_vsi_170904(parms, timing_parms, dofigs,doSub)
%%

n=1;
Ncombs = 20;
sigs = zeros(Ncombs,Ncombs);
snrs = zeros(Ncombs,Ncombs);

dofigs =0;

for TRs=linspace(2 ,5,Ncombs)
    m=1;
    for dels=linspace(0.5,2, Ncombs)
        timing_parms.t_delay = dels*ones(Nframes,1);
        timing_parms.t_adjust = TRs -timing_parms.t_aq - timing_parms.t_delay;
        
        timing_parms.t_delay = 0.1*ones(Nframes,1);
        timing_parms.t_tags = dels*ones(Nframes,1);
        timing_parms.t_adjust = TRs -timing_parms.t_aq - timing_parms.t_delay - timing_parms.t_tags;
        
        if 0.1+dels+aqwindow < TRs
            obs = gen_signals_vsi_170904(parms, timing_parms, dofigs,doSub);
            sigs(m,n) = obs(end);
            snrs(m,n) = obs(end)/sqrt(TRs);
        else
            sigs(m,n) = nan;
            snrs(m,n) = nan;
        end
  
        m = m+1;
       
    end
     n = n+1;
end

figure

subplot(211)
imagesc(sigs * 100)
colormap hot
colorbar
xlabel('TR (seconds)')
set(gca,'xtick',linspace(1,Ncombs,5))
set(gca,'xtickLabel',linspace(2, 5, 5 ) )

ylabel('PID (seconds)')
set(gca,'ytick', linspace(1,Ncombs,5))
set(gca,'ytickLabel',linspace(0.5, 2, 5 ) )
title('ASL Signal intensity')


subplot(212)
imagesc(snrs / max(snrs(:)))
colormap hot
colorbar
xlabel('TR (seconds)')
set(gca,'xtick',linspace(1,Ncombs,5))
set(gca,'xtickLabel',linspace(2, 5, 5 ) )

ylabel('PID (seconds)')
set(gca,'ytick', linspace(1,Ncombs,5))
set(gca,'ytickLabel',linspace(0.5, 2, 5 ) )
title('ASL SNR')

print -dpng SNR_timing_analysis_noArt


