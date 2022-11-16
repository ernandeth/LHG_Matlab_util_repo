%simulates an analysis a bunch of times and finds the time where the 
% maximum t score occurred.  It computes the mean and the standard 
% deviation of the peak time.

clear
tau = 6;
iti = 16;
shift = [-2:0.2:2];
avgs = 10000;
result=[];


for noise_var=0.05:0.5:4.05
    
    
    tscores=[];
    maxt=[];
    
    for i=1:avgs
        noise = randn(600,1)* sqrt(noise_var);
        %noise = randn(1200,1);
        
        t=bold_shift_sim1b( noise, tau, iti);
        tscores = [tscores; t];
        [m i]= max(t);
        maxt = [maxt; [m i]];
        fprintf('\rsimulation \t%d  var:  \t%5f',i,noise_var);
    end
    
    if noise_var==3.05
        m = maxt(:,2)*0.2 - 2.2
        hist(m,21)
        str = sprintf('Histogram of Shift Estimates (noise = %.2f)',noise_var);
        title(str,'fontsize',16);
	xlabel('Shift estimate (seconds)','fontsize',16)
	axis([-2 2 0 1000])
    end
    
    meanmaxt=mean(maxt);
    stdmaxt=std(maxt);
%     hold on
%     plot(shift, mean(tscores), 'r')
%     
%     xlabel('Time Shift (seconds)');
%     ylabel('t score');
%     str=sprintf('max t score %f +/- %f  occurs at %f +/- %f', meanmaxt(1),stdmaxt(1),meanmaxt(2)*0.1 -1.1,stdmaxt(2)*0.1);
%     
%     title(str)
    
    result = [result; noise_var  meanmaxt(1)  stdmaxt(2)*0.1];   
end

result;


figure

plot(result(:,1), result(:,3),'k','linewidth', 3)
xlabel('Noise level (var. of noise relative to BOLD)','fontsize',16)
ylabel('Std. Dev. of estimate (seconds)','fontsize',16)


save variances result m
