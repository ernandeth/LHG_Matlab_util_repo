clear
times = 1000
avg_t= zeros(21,21);
noise_var = 1.5;

all_t = [];

for count = 1:times
   
   noise = randn(600,1)* sqrt(noise_var);
   tmp = [];
   
   for i= -10:10
      t = bold_shift_sim(noise,6, 16, i*0.05);
      %t = bold_shift_data('data.dat','onsets.dat',6, i*0.05);
      tmp = [tmp ; t];
      fprintf('\rcount: %d   i: %d',count, i);  
   
   end
   
   avg_t = avg_t + tmp;
   all_t = [all_t ; tmp];
   
   
end

avg_t= avg_t /times;
description = 'Correlation as a function of time shift (-1:1) and Error in Tau (-0.25:0.25)'
save tau_error_noisy.mat avg_t all_t times description noise
%save tau_error_data.mat avg_t all_t times description noise
return

figure
t = [-1:0.1:1];
plot(t,avg_t(1,:),'--k')
hold on
plot(t,avg_t(11,:),'k')
plot(t,avg_t(21,:),'.-k')
title('Effect of Errors in the Model of the HRF the Temporal Sensitivity')
xlabel('Temporal Shift of the Regressor (seconds)')
ylabel('t score')
legend('Error in Tau = -0.5', 'No Error', 'Error in Tau = 0.5')


figure
dt = [-0.5:0.05:0.5];
plot(dt, avg_t(:,11),'k')
title('Effect of Errors in the Model of the HRF on the t score')
xlabel('Error in the Tau constants (seconds)')
ylabel('t score')

%%%%%%%%%%%%%%%%%%%%%%%%%


%x = zeros(100,4);
%x(20,1) = 1;
%x(40,2) = 1;
%x(80,4) = 1;
%y = make_hrf(1.2, 2.5, 100)
%c = conv2(x,y)
%c = c(1:100,:)


