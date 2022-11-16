function H = HRF_mat(tau1, tau2, T)
% function H = HRF_mat(tau1, tau2, T)
%
% tau1 = dispersion of the reponse (in samples),
% tau2 = dispersion of undershoot (in samples),
% T = number of points in the timecourse
%
% this function generates an impulse reponse function and uses
% it to build a system matrix so that you can quickly make HRF
% timecourses instead of convolving.  You can use that system matrix
% also to solve the minimization problems



h = inline(' (t).^(tau1).*exp(-t)/gamma(tau1) - 0.1*(t).^(tau2).*exp(-t)/gamma(tau2)',...
   'tau1','tau2','t');

% h = inline(' (t).^(tau1).*exp(-t)/gamma(tau1) - 0.3*(t).^(tau2).*exp(-t)/gamma(tau2)',...
%   'tau1','tau2','t');


t = [0:2:2*T-2]';


hh = h(tau1,tau2, t+7);

%hh = spm_hrf(1,[6 16 1 1 6 0 T-1]);

hh = hh/max(hh);
hh = hh - mean(hh);
%hh = hh/norm(hh);

% hh1 = hh;
% hh = zeros(size(t));
% hrf = spm_hrf(1);
% hh(1:length(hrf)) = hrf;

%plot(hh); hold on, plot(hh1,'r'), axis([0 30 -0.5 1])

H = toeplitz( [hh(1); zeros(T-1,1)], hh')';

H = H(:,1:T);
% for i=1:T,
%     H(:,i)=H(:,i)/norm(H(:,i));
% end
% %H = H/max(H(:));
H = H-ones(T,1)*mean(H);


return

