clear all
% low frequency, low power simulation
hh_sim_scr
count=1;
start=0;
N_cycles=25;
for Amp=linspace(0,10,25)
    for w=linspace(0,1,30)
        Iext=pulse(0:dt:T,start,N_cycles,-Amp,w);
        t=0; hh_sim_scr; 
        Vmax(count)=max(abs(V)); count = count+1;
        %drawnow; pause(0.01)
    end
end

Vmax = reshape(Vmax, 30,25);
figure; imagesc(Vmax)
ylabel('frequency (0 - 1 KHz.)'); 
xlabel('Current Amplitude (0-10 uA) ');
title('Action potentials as a function of the current amplitude and sinusoidal frequency')
print -djpeg amplitude_freq

%%

clear all
hh_sim_scr
count=1;
start=0;
w=0;
for Amp=linspace(0,10,25)
    for N_cycles = linspace(0,5,30)
        
        Iext=pulse(0:dt:T,start,N_cycles,-Amp,w);
        t=0; hh_sim_scr; 
        Vmax(count)=max(abs(V)); count = count+1;
        % drawnow; pause(0.1)
        pause
end

Vmax = reshape(Vmax, 30,25);
figure; imagesc(Vmax)
ylabel('N_cycles (0 - 5 ms)'); 
xlabel('Current Amplitude (0-10 uA) ');
title('Action potentials produced by a single pulse')
print -djpeg amplitude_N_cycles_single

%% now test N_cycless


clear all
hh_sim_scr

start=0;
w=0.3;
plotn=1;
for w=0:0.1:0.5
    count=1;
    for Amp=linspace(0,10,25)
        for N_cycles = linspace(0,10,30)
            figure(1)

            Iext=pulse(0:dt:T,start,N_cycles,-Amp,w);
            t=0; hh_sim_scr;
            Vmax(count)=max(abs(V)); count = count+1;
            % drawnow; pause(0.1)
        end
    end


    figure(33)
    subplot(2,3,plotn)
    plotn = plotn+1;
    Vmax = reshape(Vmax, 30,25);
    imagesc(Vmax)
    ylabel('N_cycles (0 - 10 ms)');
    xlabel('Current Amplitude (0-10 uA) ');
    title([ num2str(w) ' KHz sinusoid'])

end

print -djpeg amplitude_N_cycles_waves

% simulate a field from a Magstim pulse
clear all
hh_sim_scr

start=0;
N_cycles=.3
w=3.3;

count=1;
for Amp=linspace(0,1000,25)

    figure(1)

    Iext=pulse(0:dt:T,start,N_cycles,-Amp,w);
    t=0; hh_sim_scr;
    Vmax(count)=max(abs(V)); count = count+1;
    % drawnow; pause(0.1)
end
figure
Amp=linspace(0,1000,25)
plot(Amp, Vmax); title('Magstim stimulus, with increasing current')
xlabel('current [uA]'); ylabel('max Voltage across membrane')
print -djpeg magstim


%%
% higher frequency and power simulation
%
clear all
hh_sim_scr
Vss = V(end);
count=1;
start=0;
duration=5;

allAmps=linspace(0,100,25)
allFreqs=linspace(0.1,3,15)
        
for Amp=allAmps
    for w=allFreqs
        
        N_cycles = round(duration * w);
        N_cycles = 10;
        
        Iext=pulse(0:dt:T,start,N_cycles,-Amp,w);
        t=0; 
        V(:) = Vss;
        hh_sim_scr; 
        Vmax(count)= max(abs(V)); 
        NRG(count) = sum((Iext.^2)*dt);
        count = count+1;
        drawnow; pause(0.01)
    end
end




Vmax = reshape(Vmax, length(allFreqs),length(allAmps));
NRG = reshape(NRG, length(allFreqs),length(allAmps));

minNRG = zeros(size(allFreqs));
minAmp = zeros(size(allFreqs));
for count=1:length(allFreqs);
    tmp = Vmax(count,:);
    iVmax = min(find( tmp > 50));
    if ~isempty(iVmax)
        minNRG(count) = NRG(count, iVmax);
        minAmp(count) = allAmps(iVmax);
    end
end

figure; imagesc(Vmax)
ylabel('frequency (0 - 10 KHz.)'); 
xlabel('Current Amplitude (0-100 uA) ');
title('Max. Voltage as a function of the current amplitude and sinusoidal frequency')


figure(33)
subplot(211)
plot(allFreqs, minAmp) ; title('minimum Current Amplitude required per frequency')
ylabel('I_{min} (uA)')
subplot(212)
plot(allFreqs, minNRG) ; title('minimum ENERGY required per frequency'); 
xlabel('Frequency (kHz)')
ylabel('E_{min} (uA^2*ms)')




%%
% higher frequency and power simulation
% as a function of the number of cycles
clear all
hh_sim_scr
count=1;
start=0;
duration=5;

allAmps=linspace(0,100,25)
allFreqs=linspace(0,3,30)
allN = 1:30;

w = 1;

for Amp=allAmps

    for N_cycles= allN
        
        Iext=pulse(0:dt:T,start,N_cycles,-Amp,w);
        t=0; hh_sim_scr; 
        Vmax(count)= max(abs(V)); 
        NRG(count) = sum((Iext.^2)*dt);
        count = count+1;
        drawnow; pause(0.01)
    end
end

Vmax = reshape(Vmax, length(allN),length(allAmps));
NRG = reshape(NRG, length(allN),length(allAmps));

minNRG = zeros(size(allN));
minAmp = zeros(size(allN));
for count=1:length(allN);
    tmp = Vmax(count,:);
    iVmax = min(find( tmp > 50));
    if ~isempty(iVmax)
        minNRG(count) = NRG(count, iVmax);
        minAmp(count) = allAmps(iVmax);
    end
end

figure; imagesc(Vmax)
ylabel('Number of cycles '); 
xlabel('Current Amplitude (0-100 uA) ');
title('Maximum Voltage as a function of the current amplitude and N. cycles')
print -djpeg amplitude_1KHz

figure
subplot(211)
plot(allN, minAmp) ; title('Minimum Current Amplitude required ')
ylabel('I_{min} (uA)')
subplot(212)
plot(allN, minNRG) ; title('Minimum ENERGY required '); 
xlabel('Number of cycles Used')
ylabel('E_{min} (uA^2*ms)')



