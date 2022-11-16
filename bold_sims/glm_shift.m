clear
times = 10;

avg_t_glm = zeros(5,21);

noise_level = 3;
NSAMPLES = 1200;

all_t = [];
tic
for i=1:times
   toc
   %noise = noise_level *randn(NSAMPLES,1);   
   noise = noise_level *make1overf(NSAMPLES,0.5);   
   beta = abs(6* randn(5, 1)) + 1;
   
   %bold_shift_sim(noise, 1.5, 4, -0.2)
   
   t_score = fake_spm_taus(-0.2, noise, beta);
   avg_t_glm(1,:) = avg_t_glm(1,:) + t_score;
   all_t = [all_t ; t_score];
   
   
   
   t_score = fake_spm_taus(-0.1, noise, beta);
   avg_t_glm(2,:) = avg_t_glm(2,:) + t_score;
   all_t = [all_t ; t_score];
   
   
   
   
   t_score = fake_spm_taus(0, noise, beta);
   avg_t_glm(3,:) = avg_t_glm(3,:) + t_score;
   all_t = [all_t ; t_score];
   
   
   
   t_score = fake_spm_taus(0.1, noise, beta);
   avg_t_glm(4,:) = avg_t_glm(4,:) + t_score;
   all_t = [all_t ; t_score];
   
   
   
   t_score = fake_spm_taus(0.2, noise, beta);
   avg_t_glm(5,:) = avg_t_glm(5,:) + t_score;
   all_t = [all_t ; t_score];
   
   i
   
   %noise_level = noise_level + 3   
   
end   

avg_t_glm = avg_t_glm/times

description = 'GLM as a function of time shift(-1:1). 5 regressors. SPM''s HRF. Tau_error=0.1'

save glm_shift_error2_noisefree

shift = [-1:0.1:1];

plot(shift, avg_t_glm)
hold on
plot(shift, avg_t_glm(3,:), 'r')

xlabel('Time')
ylabel('t score')
grid

title('Sensitivity of the GLM to time shifts of one of the regressors when the parameters of HRF are wrong. SPM hrf')

hold off