clear
times = 1000 
avg_t= zeros(11,21);
noise_var = 1.5;

all_t = [];

for count = 1:times
   
   noise = randn(600,1)* sqrt(noise_var);
   tmp = [];
   
   for i= 0:0.1:1
      t = bold_shift_sim2(noise,6, 16, i);
      %t = bold_shift_data2('data.dat','onsets.dat',6,i);

      tmp = [tmp ; t];
      fprintf('\rcount: %d   i: %d',count, i);  
   
   end
   
   avg_t = avg_t + tmp;
   all_t = [all_t ; tmp];
   
   
end

avg_t= avg_t /times
description = 'Correlation as a function of time shift (-1:1) and Error in ITI (-1:1)'
save iti_error_noisy.mat avg_t all_t times description noise
%save iti_error_data.mat avg_t all_t times description noise
%save iti_error.mat avg_t all_t times description noise

figure
axes('fontsize',8);
t = [-1:0.1:1];
plot(t,avg_t(1,:),'k')
hold on
plot(t,avg_t(11,:),'--k')

title('Effect of Errors in the Timing of Events on the Temporal Sensitivity', 'fontsize',8)
xlabel('Temporal Shift of the Regressor (seconds)', 'fontsize',8)
ylabel('t score', 'fontsize',8)
legend('No ITI Error', 'ITI error jitter = 1')

figure
axes('fontsize',8);
d_iti = [0:0.1:1];
plot(d_iti, avg_t(:,11),'k')
title('Effect of Errors in the Timing of Events on the t score', 'fontsize',8)
xlabel('Variance of the ITI Error (seconds)', 'fontsize',8)
ylabel('t score', 'fontsize',8)
