clear
times = 1000
avg_t= zeros(21,21);
noise_var = 0.001;
all_t = [];


for count = 1:times
   
   noise = randn(600,1)* sqrt(noise_var);
   tmp = [];
   
   for i= -10:1:10
      fprintf('count: %d   i: %d  ',count, i);  
      t = bold_shift_sim3(noise,6, 16, i);
      %t = bold_shift_data3('data.dat','onsets.dat', 6, i);
      tmp = [tmp ; t];
   end
      
   avg_t = avg_t + tmp;
   %all_t = [all_t ; tmp];
   
   
end

avg_t= avg_t /times
description = 'Correlation as a function of time shift (-1:1) and Error in NUmber of events (-10:10)'
%save n_iti_error_noisy.mat avg_t all_t times description noise noise_var
save n_iti_error.mat avg_t all_t times description noise noise_var
%save n_iti_error_data.mat avg_t all_t times description noise noise_var

figure
t = [-1:0.1:1];
plot(t,avg_t(1,:),'-.k')
hold on
plot(t,avg_t(11,:),'k')
plot(t,avg_t(21,:),'--k')
title('Effect of Errors in the Number of Trials on the Temporal Sensitivity')
xlabel('Temporal Sensitivity (seconds)')
ylabel('t score')
legend('Number of Extra Events = -10', ...
   'perfect match', ...
   'Number of extra events = 10')

figure
d_niti = [-10:10];
plot(d_niti, avg_t(:,11),'k')
title('Effect of Errors in the Number of Events on the t score' )
xlabel('Number of Errors')
ylabel('t score')
