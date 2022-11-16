%%%%%%%%%%%%%%
% compute the numerical simulation using the multivariate regression approach
% Compute the analytical simulation of the time shifts
% do both for different amounts of noise
clear
tau = 6;
iti = 16;
shift = [-1:0.1:1];
avgs = 1000;
df = 1198;
noise_var=0.5;


%%%%%%%%%%%%%%
% compute the numerical simulation using the multivariate regression approach
%figure
axes('fontsize',8);

corr_glm = zeros(10,21);

for i=1:avgs
   noise = randn(600,1);

   c=fake_spm( noise, tau, iti, 1, [1]);
   
   corr_glm = corr_glm + c;
   i
end

corr_glm = corr_glm/avgs;

subplot 121, plot(shift, corr_glm, 'k')
xlabel('Time Shift (seconds)' );
ylabel('t score' );
title('Monte Carlo Simulation' )

%%%%%%%%%%%%%%
% Compute the analytical simulation of the time shifts


%figure
axes('fontsize',8);
noise_var = 0.5;

corr_ana = bold_shift_ana10(noise_var,tau,iti);
corr_ana = corr_ana';

%t_score = corr_ana * sqrt(df) ./ sqrt(1 - corr_ana.^2)
%corr_ana = t_score

subplot 122, plot(shift, corr_ana, 'k')
xlabel('Time Shift (seconds)');
ylabel('t score' );
title('Analytical expression' )


return
%%%%%%%%%%
% Compute the numerical simulation using a correlation coefficient.


corr_sim = zeros(10,21);
%noise_var = 3;

for i=1:avgs
   %noise = randn(1200,1)* sqrt(noise_var);
   noise = randn(600,1);
   c=bold_shift_sim(noise, tau, iti, 0);
   corr_sim = corr_sim + c;
   i
end

corr_sim = corr_sim/avgs;


%t_score = corr_sim * sqrt(df) ./ sqrt(1 - corr_sim.^2)
%corr_sim = t_score
%figure
axes('fontsize',8);

subplot 131, plot(shift, corr_sim, 'k')
xlabel('Time Shift (seconds)' );
ylabel('t score' );
title('Correlation' )
%%%%%%%%%%%%%%%%%%%%%%%%

save simulate.mat corr_sim corr_ana corr_glm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the effect of time shifts with the actual Garavan experiment 
% event times read from a file.  Using the GLM.

figure
corr_glm_garavan = zeros(1,21);

for i=1:avgs
   noise = randn(1200,1)* noise_var+100000;
   %noise = randn(1200,1)* sqrt(noise_var)+100000;
   %noise = make1overf(1200, 1) * sqrt(noise_var);
   c=fake_garavan( noise, tau, iti, 1, [0 0 0 1 -1 ]);
   corr_glm_garavan = corr_glm_garavan + c;
   i
end

corr_glm_garavan = corr_glm_garavan/avgs;

plot(shift, corr_glm_garavan, 'k')
xlabel('Time Shift (seconds)');
ylabel('t score');
title('Simulated Experiment with five regressors')
