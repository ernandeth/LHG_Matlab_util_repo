vel = [1 5 10 15 20 30]
alphaNoMT = zeros(1,5);
alphaMT = zeros(1,5);
alphaNoMTend = zeros(1,5);
alphaMTend = zeros(1,5);
alpha_relaxation_ratio = zeros(1,5);
alpha_reduction = zeros(1,5);

for v=1:length(vel)
	mt = load(sprintf('MTVel%02d.dat', vel(v)));
	nomt = load(sprintf('NOMTVel%02d.dat', vel(v)));
	subplot (121), plot(nomt(:,1), nomt(:,2)), hold on
	axis([-1 1 -1 1]);

        string1 = sprintf('MTVel%02d.dat', v);
	hold on; grid on;
        subplot 122, plot(mt(:,1),mt(:,2),'k'), hold on
        title(string1);
	axis([-1 1 -1 1]);
	drawnow

	alphaNoMT(v) = min(nomt(:,2));
 	alphaNoMT(v) = (1 - alphaNoMT(v))/2;
   
 	alphaMT(v) = min(mt(:,2));
 	alphaMT(v) = (1 - alphaMT(v))/2;
    
	alphaNoMTend(v) = (1-nomt(end,2))/2;
	alphaMTend(v) = (1-mt(end,2))/2;
    
	alpha_reduction(v) = 100*(alphaNoMT(v) - alphaMT(v))/alphaNoMT(v) 
	alpha_relaxation_ratio(v) = 100*((alphaNoMT(v)-alphaNoMTend(v)) - (alphaMT(v)-alphaMTend(v))) / (alphaNoMT(v)-alphaNoMTend(v))

end
title('B. Effect of Velocity (MT present)');
xlabel('Time (s)'), ylabel('Mz/Mz0')
dofontsize(16)
fatlines
title('B. Effect of Velocity (MT present)');
xlabel('Time (s)'), ylabel('Mz/Mz0')
dofontsize(16)
fatlines

mt = load ('MTVel20.dat');
subplot 121, hold off
plot(mt(:,1),mt(:,2),'k')
title('A. Effect of Magnetization Transfer');
    
nomt = load ('NOMTVel20.dat');
hold on;     grid on;
plot(nomt(:,1),nomt(:,2),'--k')
axis([-1 1 -1 1]);
legend('MT present' ,'No MT')
xlabel('Time(s.)')
ylabel('Mz (a.u)')
fatlines
dofontsize(16)
    
alpha_relaxation_ratio = ...
    	(alphaMT - alphaMTend)./alphaMT*100  - ...
	(alphaNoMT - alphaNoMTend)./alphaNoMT*100  

% net loss of alpha due to MT:
100*(alphaMTend - alphaNoMTend )./ alphaNoMTend

figure
plot(vel, alpha_reduction)	
figure
plot(vel, alpha_relaxation_ratio)	
