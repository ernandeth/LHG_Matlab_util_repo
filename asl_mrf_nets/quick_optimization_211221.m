
dim1 = 1;
dim2 = 1;
DEBUG = 1;

Ntests = 8;
all_NRMS = zeros(Ntests);
all_SNR_eff = zeros(Ntests);
all_seqDuration = zeros(Ntests);
all_Ncycles = zeros(Ntests);
all_maxtimes = zeros(Ntests);

maxtime = 1.6;
Ncycles = 15;
best_SNR_eff = 0;
Nframes = 250;
Nbreak = 7;

Ncycle_vals = linspace(1,20, Ntests);
maxTime_vals =  linspace(0.1, 1.4, Ntests);

% after optimizing:
Ncycle_vals = 5
maxTime_vals = 0.5

parms.mtis0 =     1 ;
parms.Disp =      40;
parms.r1blood = 1/1.7;

parms.f =       0.01;
parms.cbva =    0.01;
parms.bat =     0.5;
parms.r1tis =   1/0.9;%1.4;
parms.flip =    deg2rad(40);
parms.r2tis=    1/0.090;


parms2 = parms;
 
%parms2.f = 0.003;
parms2.bat = 0.3;
parms2.cbva = 0.02;
%parms2.r2tis = 1/0.05;

for Ncycles = Ncycle_vals
    dim1 = 1;
    for maxtime = maxTime_vals
        
        
        
        %figures and subtraction
        doFigs = 0;
        doSub = 0;
        
        best_NRMS=0;
        
        % timing parms
        %timing_parms.label_type = 'FTVSI-sinc';%'FTVSI-sinc'; % 'BIR8inv'; % 'BIR8'
        timing_parms.label_type = 'BIR8'; % 'BIR8'
        %timing_parms.label_type = 'BIR8inv'; % 'BIR8'
        timing_parms.readout_type = 'FSE'; %'GRE';
        
        %timing_parms.t_delay = linspace(0.2,1.6, Nframes)';
        %timing_parms.t_delay                = [linspace(0.2, 0.3, floor(Nframes/4)) ...
        %    linspace(1, 1.4, Nframes/2) ...
        %    linspace(0.3, 0.4, ceil(Nframes/4))]';
        
        
        
        % timing_parms.t_delay = 0.2 + (maxtime/2)*(1+chirp(linspace(1,0,Nframes), 0, 1, Ncycles))';
        % timing_parms.t_adjusts = 0.2 + (maxtime/2)*(1+chirp(linspace( 0,1, Nframes), 0, 1,Ncycles))';
        %
        %     timing_parms.t_delay(:) = maxtime;
        %     timing_parms.t_adjusts(:) = 3-maxtime;
        
        %{
        % Anish schedule:
        tmp = load('timing_files_201/t_adjusts.txt');
        timing_parms.t_adjusts = tmp(:);
        tmp = load('timing_files_201/t_delays.txt');
        timing_parms.t_delay = tmp(:);
        %}
        
        
        timing_parms.t_tag                  = 0.050 *ones(Nframes,1) ;
        timing_parms.ArtSup_delay           = 0.150 *ones(Nframes,1) ;
        timing_parms.t_aq                   = 0.0329 *18*ones(Nframes,1) ;
        timing_parms.t_aq                   = 0.42938*ones(Nframes,1) ;
        timing_parms.order                  = 1;
        
        % LHG 3/15/22
        timing_parms.t_delay                = 0.25 + (maxtime/2)*(1+chirp(linspace(0,1, Nframes), 0, 1, Ncycles))';
        timing_parms.t_adjusts              = 0.5*ones(Nframes,1) ;
        timing_parms.ArtSup_delay           = 0.05 + (maxtime)*(1+chirp(linspace( 1,0, Nframes), 0, 1,1.25*Ncycles))';
        
        delays = timing_parms.t_delay;
        ordr = timing_parms.order;
        
        % Make timing paramaters:
        best_NRMS = 0;
        for t=1:1  % random orders : choose the best_NRMS
            
            % random distribution of control and label pulses
            rp  = randperm(Nframes);
            rp=rp(1:end/2);
            % random distribution of arterial suppression (1/4 of the time)
            rp2 = randperm(Nframes);
            rp2=rp2(1:end/2);
            % random distribution of delays:
            rp3 = randperm(Nframes);
            %
            timing_parms.labelcontrol       = ones(Nframes,1);
            timing_parms.labelcontrol(rp)   = 0;
            
            timing_parms.doArtSup           = ones(Nframes,1);
            timing_parms.doArtSup(rp2)      = 0;

            timing_parms.doArtSup(1:Nbreak:end)      = -1;
            timing_parms.labelcontrol(5:Nbreak:end)   = -1 ;
         
         % generate Data : single run
            test_signal = abs(single(  ...
                gen_signals_vs_220421( parms,...
                timing_parms,...
                doFigs,doSub)));% ...
            %+ (0*(randn(1,Nframes))))); % LH Model one BAT
            
            
            test_signal2 = abs(single(  ...
                gen_signals_vs_220421( parms2,...
                timing_parms,...
                doFigs,doSub)));% ...
            
            rms_change = norm(test_signal-test_signal2)/norm(test_signal);
            if rms_change > best_NRMS
                rp_best= rp;
                rp2_best = rp2;
                best_NRMS = rms_change;
            end
            
        end
        
        % now test the signal wit the best RMS
        %
        timing_parms.labelcontrol           = ones(Nframes,1);
        timing_parms.labelcontrol(rp_best)  = 0;
        timing_parms.doArtSup                = ones(Nframes,1);
        timing_parms.doArtSup(rp2_best)      = 0;

         timing_parms.doArtSup(1:Nbreak:end)      = -1;
         timing_parms.labelcontrol(5:Nbreak:end)   = -1 ;
%        timing_parms.labelcontrol(:)   = -1 ;

        %}
        
        seqDuration = sum(timing_parms.t_delay)+ ...
            sum(timing_parms.t_tag) + ...
            sum(timing_parms.t_adjusts) + ...
            sum(timing_parms.ArtSup_delay)+...
            sum(timing_parms.t_aq)
        
        SNR_efficiency = 1e4*best_NRMS/sqrt(seqDuration);
        
        if (best_SNR_eff < SNR_efficiency)
            best_SNR_eff = SNR_efficiency;
            best_timing_parms = timing_parms;
        end
        
        all_Ncycles(dim1, dim2) = Ncycles;
        all_maxtimes(dim1, dim2) = maxtime;
        
        all_NRMS(dim1, dim2) = best_NRMS;
        all_SNR_eff(dim1, dim2) = SNR_efficiency;
        all_seqDuration(dim1, dim2) = seqDuration;
        
        fprintf('\nmaxTime: %0.2f sec., \tNcycles: %0.2f, \tNRMS: %0.2f , \tSNR eff: %0.2f ',...
            maxtime, Ncycles, best_NRMS, SNR_efficiency)
        
        dim1 = dim1+1;
        
        % generate Data : single run
        %parms.flip = deg2rad(90);
        
        test_signal = abs(single(  ...
            gen_signals_vs_220421(parms,...
            timing_parms,...
            doFigs,doSub))); % ...
        %+ (0*(randn(1,Nframes))))); % LH Model one BAT
        
        test_signal2 = abs(single(  ...
            gen_signals_vs_220421(parms2,...
            timing_parms,...
            doFigs,doSub))); ...
            %+ (0*(randn(1,Nframes))))); % LH Model one BAT
        %
        if DEBUG
            figure(1)
            subplot(211)
            plot(timing_parms.t_delay); title('t delay')
            axis tight
            
            subplot(212)
            plot(timing_parms.ArtSup_delay); title('AS delay')
            axis tight
            
            figure(2)
            
            subplot(211)
            plot(test_signal); hold on;
            plot(test_signal2); hold off
            
            subplot(212)
            plot( (test_signal-test_signal2)/norm(test_signal) )
            title(sprintf('NRMS change  : %0.2e',...
                best_NRMS));
            %            axis([1 length(test_signal) 0 1e-3])
            
            text(10, double(min(test_signal2)), sprintf('Duration= %0.2f seconds', seqDuration));
            drawnow
        end
        
        
        %}
    end
    
    dim2 = dim2+1;
end
figure
intervals = timing_parms.t_delay+timing_parms.ArtSup_delay;
plot(intervals)
title('AStime + t delay')

%make_schedule_files(timing_parms,'chirp_vssx2_schedule_250' )
%save optimal_sequence.mat best_timing_parms


%%
%{
figure(3)
subplot(223)
imagesc(all_SNR_eff)
xlabel('N cycles')
ylabel('Max. Time')
title('SNR efficiency for \Delta CBV')
colormap gray
colorbar

xticklabels(string(round(100*Ncycle_vals)/100))
yticklabels(string(round(100*maxTime_vals)/100))
drawnow

subplot(221)
imagesc(all_NRMS)
xlabel('N cycles')
ylabel('Max. Time')
title('NRMS for \Delta CBV')
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
