% simulate a time series of VS ASL
% with different PLD's and TRs
% compute the SNR per unit time
close all

TR = 4;
doSub = 0;
dofigs = 1;


parms.f=          0.01 ;
parms.cbva =      0.02;
parms.bat =       0.13 ;
parms.r1tis =     1/1.4  ;
parms.Mtis0 =     1 ;
parms.flip =      90*pi/180 ; % flip angle in radians
parms.r2tis =     1/0.090 ;

% include a B1 error term in the labeling pulses
% 0 error means multiplying by 1.
% 1% error mean multiplying by 1.01
parms.b1err =     0.05;


% pulse sequence parameters
aqparms.Nframes = 6;
aqparms.label_type =    'BIR8inv'; %'FTVSS'; %'FTVSI-sinc'; % 'BIR8inv'; % 'BIR8'
aqparms.RO_type =       'FSE'; %'FSE';   % 'GRE'
aqparms.t_tags =        0;% 0.1*ones(Nframes,1);
aqparms.del2 =          1.8*ones(aqparms.Nframes,1);
aqparms.del3 =          .15*ones(aqparms.Nframes,1);  % delay between AS pulse and acqusition
aqparms.labelcontrol =  zeros(aqparms.Nframes,1);
aqparms.labelcontrol(2:2:end)= 1;
aqparms.labelcontrol(1:2:end) =   0;
aqparms.order =         1;
aqparms.doArtSup =      ones(aqparms.Nframes,1);
aqparms.t_aq =          0.60 ;  % duration of the whole readout
aqparms.Nkz =           16;
aqparms.del1 =          TR - aqparms.del2 -  aqparms.del3 -  aqparms.t_aq;

obs = gen_signals_vs_230718(parms, aqparms, dofigs,doSub)
%%
n=1;
Ncombs = 50;
sigs = zeros(Ncombs,Ncombs);
snrs = zeros(Ncombs,Ncombs);

dofigs =0;

dels=linspace(0.1,2, Ncombs)
TRs=linspace(2 ,5,Ncombs)

for tmpTR = TRs
    m=1;
    for tmpDel=dels
        aqparms.del2 = tmpDel*ones(aqparms.Nframes,1);
        if aqparms.t_aq + aqparms.del2 + aqparms.del3 < tmpTR
            aqparms.del1 = tmpTR -aqparms.t_aq - aqparms.del2-aqparms.del3;
            obs = gen_signals_vs_230718(parms, aqparms, dofigs,doSub, 1e-3);
            obs = abs(obs(1:2:end)-obs(2:2:end));
            sigs(m,n) = obs(end);
            snrs(m,n) = obs(end)/sqrt(tmpTR);
            drawnow
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
set(gca,'xtickLabel',linspace(min(TRs), max(TRs), 5 )  )

ylabel('PID (seconds)')
set(gca,'ytick', linspace(1,Ncombs,5))
set(gca,'ytickLabel',linspace(min(dels), max(dels), 5 ) )
title('ASL Signal intensity')


subplot(212)
imagesc(snrs / max(snrs(:)))
colormap hot
colorbar
xlabel('TR (seconds)')
set(gca,'xtick',linspace(1,Ncombs,5))
set(gca,'xtickLabel',linspace(min(TRs), max(TRs), 5 )  )

ylabel('PID (seconds)')
set(gca,'ytick', linspace(1,Ncombs,5))
set(gca,'ytickLabel',linspace(min(dels), max(dels), 5 ) )
title('ASL SNR')

print -dpng SNR_timing_analysis_noArt


