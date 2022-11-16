N = 258;

VarR = 100;  % Ratio of WN variance to AR process variance

V=(spm_Q(0.95,N)+eye(N)*VarR)/(1+VarR);
Spec=abs(fft(full(V(:,floor(N/2))))).^2;
% Select first half - 
if rem(N,2),         % nfft odd
  select = (1:(N+1)/2+1)';  % don't include DC or Nyquist components
else
  select = (1:(N)/2+1)';
end
      
% Calculate the single-sided spectrum which includes the full power
Pwr = [2*Spec(select,:)];
Freq = (select-1)/max(select-1)*0.5;

plot(Freq,Pwr)

title(sprintf('AR+WN power    \\rho=%g    (VarWN)/(VarAR) = %g',rho,VarR))
print -dtiff ARWNpower.tif