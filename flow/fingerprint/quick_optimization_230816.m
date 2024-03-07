
dim1 = 1;
dim2 = 1;
DEBUG = 0;


maxtime = 1;
Ncycles = 5;
best_RMSdiff_eff = 0;
Nframes = 400;

Ntests = 20;
Ncycle_vals = linspace(0,16, Ntests);
maxTime_vals =  linspace(0.5, 1.0, Ntests);

all_RMS = zeros(Ntests);
all_RMSdiff_eff = zeros(Ntests);
all_seqDuration = zeros(Ntests);
all_Ncycles = zeros(Ntests);
all_maxtimes = zeros(Ntests);

% after optimizing:
% Ncycle_vals = 5
% maxTime_vals = 0.5
%
Nechoes = 1;

parms.Mtis0 =     1 ;
parms.Disp =      40;
parms.r1blood = 1/1.7;

parms.f =       0.01;
parms.cbva =    0.01;
parms.bat =     0.001;
parms.r1tis =   1/1.2;
parms.flip =    deg2rad(30);
parms.r2tis=    1/0.090;
parms.b1err = 0;


parms2 = parms;

parms2.f = 0.02;
parms2.bat = 0.2;

% periodicity for label and control
n=linspace(0,2*pi, Nframes)*20;
Nframes00 = 25;
Nbreaks = 10; % don't so any prep pulses every Nbreaks  
best_diff = 0;

for Ncycles = Ncycle_vals
    dim1 = 1;
    for maxtime = maxTime_vals
        %figures and subtraction
        doFigs = 0;
        doSub = 0;

        best_RMS=0;

        % timing parms
        %timing_parms.label_type = 'FTVSI-sinc';%'FTVSI-sinc'; % 'BIR8inv'; % 'BIR8'
        timing_parms.label_type              = 'BIR8'; % 'BIR8'
        %timing_parms.label_type = 'BIR8inv'; % 'BIR8'
        timing_parms.readout_type           = 'GRE'; %'GRE';
        timing_parms.RO_type                = 'GRE';
        timing_parms.t_tags                 = 0.050 * ones(Nframes,1) ;
        timing_parms.RO_time                = 0.0329 *Nechoes*ones(Nframes,1) ;
        timing_parms.t_aq                   = timing_parms.RO_time;
        timing_parms.order                  = 1;
        timing_parms.flip                   = deg2rad(30);

        % LHG 12/15/23
        timing_parms.del1  = 0.1 + 0*(maxtime/2)*(1+chirp(linspace(0,1, Nframes), 0, 1, Ncycles, 'quadratic'))';
        timing_parms.del2  = 0.05 +   (maxtime/2 )*(1+chirp(linspace(0,1, Nframes), 0, 1, Ncycles, 'quadratic'))';% .* linspace(1,2,Nframes)' ;
        timing_parms.del3  = 0.05 +    (maxtime/2/5 )*(1+chirp(linspace(1,0, Nframes), 0, 1, Ncycles, 'quadratic'))';


        delays = timing_parms.del2;
        ordr = timing_parms.order;

        % Loop over randomized  control-label order:
        best_RMS = 0;
        for t=1:1  % random prep selectivivity order : choose the best_RMS

            % randomize the order of selecetive and non selecitve once
            % random distribution of control and label pulses
            %rp  = randperm(Nframes);
            %rp=rp(1:end/2);
            rp = (square(n)'+1)/2;

            % random distribution of arterial suppression
            %rp2 = randperm(Nframes);
            %rp2 = rp2(1:end/2);
            rp2 = (square(n)'+1)/2;

            % random distribution of breaks for pulse 1 and pulse 2
            rp3 = randperm(Nframes);
            rp3 = rp3(1:Nbreaks);
            rp4 = randperm(Nframes);
            rp4 = rp4(1:Nbreaks);

            % ... but dont' take breaks in the first Nframes00 frames
            rp3 = rp3(find(rp3>Nframes00));
            rp4 = rp4(find(rp4>Nframes00));

            %

            timing_parms.labelcontrol       = rp;
            %timing_parms.labelcontrol      = ones(Nframes,1);
            %timing_parms.labelcontrol(rp)   = 0;
            %timing_parms.labelcontrol(1:2:end)   = 0;

            timing_parms.doArtSup           = rp2;
            %timing_parms.doArtSup           = ones(Nframes,1);
            %timing_parms.doArtSup(rp2)      = 0;
            %timing_parms.doArtSup(1:2:end)      = 0;

            %
            timing_parms.doArtSup(1:Nbreaks:end)      = -1;
            timing_parms.labelcontrol(1:Nbreaks:end)   = -1 ;
            %
            %{
            % Include the breaks here:
            timing_parms.doArtSup(rp3)      = -1;
            timing_parms.labelcontrol(rp4)   = -1 ;
            %}
            % generate Data : single run
            test_signal = abs(single(  ...
                gen_signals_vs_230918( parms,...
                timing_parms,...
                doFigs,doSub)));% ...

            test_signal2 = abs(single(  ...
                gen_signals_vs_230918( parms2,...
                timing_parms,...
                doFigs,doSub)));

            rms_change = norm(test_signal-test_signal2) / mean(test_signal);

            % choose the best of the random orders for this timing
            % combination:
            if rms_change > best_RMS
                bestRandparms = timing_parms;
                best_RMS = rms_change;
            end
        end

        % now test the signal with the best RMS of the control-label
        % ordering
        timing_parms = bestRandparms;

        seqDuration = sum(timing_parms.del1)+ ...
            sum(timing_parms.del2) + ...
            sum(timing_parms.del3)+...
            sum(timing_parms.RO_time);

        % RMSdiff_eficiency is the main measure of goodness
        RMSdiff_efficiency = best_RMS/sqrt(seqDuration);

        all_Ncycles(dim1, dim2) = Ncycles;
        all_maxtimes(dim1, dim2) = maxtime;

        all_RMS(dim1, dim2) = best_RMS;
        all_RMSdiff_eff(dim1, dim2) = RMSdiff_efficiency;
        all_seqDuration(dim1, dim2) = seqDuration;

        dim1 = dim1+1;

        % update the best set of parameters based on RMSdiff efficiency
        if (best_RMSdiff_eff < RMSdiff_efficiency)
            best_RMSdiff_eff = RMSdiff_efficiency;

            fprintf(['\nBEST RMSdiff efficiency so far:' ...
                '\nmaxTime: %0.2f sec., \nNcycles: %0.2f,'...
                '\nRMS change: %0.2e , \nRMSdiff efficiency: %0.2e '...
                '\nTotal Duration: %0.2f \n'],...
                maxtime, Ncycles, best_RMS, RMSdiff_efficiency,seqDuration);

            opt_timing_parms = timing_parms;
            best_diff = test_signal-test_signal2;
        end
        %%
        if DEBUG
            % Plot the best of the random orders
            test_signal = (single(  ...
                gen_signals_vs_230918(parms,...
                timing_parms,...
                doFigs,doSub))); % ...
            %+ (0*(randn(1,Nframes))))); % LH Model one BAT

            test_signal2 = (single(  ...
                gen_signals_vs_230918(parms2,...
                timing_parms,...
                doFigs,doSub))); ...
                %+ (0*(randn(1,Nframes))))); % LH Model one BAT
            %
            figure(1)
            subplot(311)
            plot(timing_parms.del1); title('del1')
            axis tight

            subplot(312)
            plot(timing_parms.del2); title('del2')
            axis tight

            subplot(313)
            plot(timing_parms.del3); title('del3')
            axis tight

            figure(2)

            subplot(211)
            plot(test_signal); hold on;
            plot(test_signal2); hold off

            subplot(212)
            hold off
            plot(best_diff, 'g')
            hold on
            plot( (test_signal-test_signal2) )
            title(sprintf('RMS change  : %0.2e',...
                best_RMS));
            legend('Best', 'Current')

            title(sprintf('Currently: %0.2f cycles , %0.2f s..', Ncycles, maxtime))
            drawnow

        end
        %%
    end

    dim2 = dim2+1;
end
%%
figure
intervals = opt_timing_parms.del1 + opt_timing_parms.del2 + opt_timing_parms.del3 + opt_timing_parms.RO_time;
plot(intervals)
title('TR schedule (sec)')

figure
subplot(223)
imagesc(all_RMSdiff_eff)
title('RMSdiff efficiency')
xlabel('Chirp Cycles')
ylabel('Chirp. Amplitude (sec)')
xticks(1:2:Ntests);
xticklabels(Ncycle_vals(1:2:end));
yticks(1:2:Ntests);
yticklabels(maxTime_vals(1:2:end));
colorbar

subplot(222)
imagesc(all_RMS)
title('RMS')
xlabel('Chirp Cycles')
ylabel('Chirp. Amplitude (sec)')
xticks(1:2:Ntests);
xticklabels(Ncycle_vals(1:2:end));
yticks(1:2:Ntests);
yticklabels(maxTime_vals(1:2:end));
colorbar

subplot(221)
imagesc(all_seqDuration)
title('seq. duration')
xlabel('Chirp Cycles')
xticks(1:2:Ntests);
xticklabels(Ncycle_vals(1:2:end));
yticks(1:2:Ntests);
yticklabels(maxTime_vals(1:2:end));
colorbar

ylabel('Chirp. Amplitude (sec)')


%%
figure
quick_mrf_sensitivity(opt_timing_parms, parms)

%
figure
make_schedule_files(...
    opt_timing_parms, ...
    './chirp_vssx2_schedule_400');
%

%save optimal_sequence.mat opt_timing_parms
%%
%{
figure(3)
subplot(223)
imagesc(all_RMSdiff_eff)
xlabel('N cycles')
ylabel('Max. Time')
title('RMSdiff efficiency for \Delta CBV')
colormap gray
colorbar

xticklabels(string(round(100*Ncycle_vals)/100))
yticklabels(string(round(100*maxTime_vals)/100))
drawnow

subplot(221)
imagesc(all_RMS)
xlabel('N cycles')
ylabel('Max. Time')
title('RMS for \Delta CBV')
colormap gray
colorbar
xticklabels(string(round(100*Ncycle_vals)/100))
yticklabels(string(round(100*maxTime_vals)/100))

subplot(222)
imagesc(all_seqDuration)
xlabel('N cycles')
ylabel('Max. Time')
title('SeqDuration for \Delta CBV')

xticklabels(string(round(100*Ncycle_vals)/100))
yticklabels(string(round(100*maxTime_vals)/100))
colorbar
drawnow
%}
