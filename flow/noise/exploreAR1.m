function exporeAR1
close all
rho=0.96;
TR=1.4;
nyq=1/(2*TR);

load ('voxels.mat');
tlen = size(mytimeseries_act,1);
NITER = size(mytimeseries_act,2);
w = (-pi: 2*pi/tlen:pi);
w = w(1:end-1);
act_data = mean(mytimeseries_act,2);
rest_data = mean(mytimeseries_rest,2);
% 
% rest_data = mytimeseries_rest(:,10);
% act_data = mytimeseries_act(:,55);

act_fdata = abs(fftshift(fft(act_data)));
rest_fdata = abs(fftshift(fft(rest_data)));
% act_fdata = circshift(act_fdata,-1);
% rest_fdata = circshift(rest_fdata,-1);

plot(w, act_fdata,'r');hold on
plot(w, rest_fdata,'g');
axis([-4 4 0 3500])
title('spectrum of the mean')
xlabel('frequency')
legend('active', 'resting')
hold off

% rest_ss = pairwise(rest_data);
% act_ss = pairwise(act_data);
% act_fss = abs(fftshift(fft(act_ss)));
% rest_fss = abs(fftshift(fft(rest_ss)));
% 
% w = (-pi/2: 2*pi/tlen:pi/2);
% w = w(1:end-1);
% 
% figure
% plot(w, act_fss,'r');hold on
% plot(w, rest_fss,'g');
% axis([-4 4 0 1500])
% title('spectrum of the pairwise subtraction')
% xlabel('frequency')
% legend('active', 'resting')
% hold off

% remove the whiteness?
rest_fdata = rest_fdata - mean(rest_fdata(10:40));

guess0 = [1e4 0.9 ];
LB = [0 0 ]; UB=[max(rest_fdata), 1];consts=[];
optvar=optimset('lsqnonlin');

% remove perfusion baseline 
data(1)=0;

%w = (linspace(-pi,pi,tlen))';
guess = lsqnonlin(@AR1spec_lsq, guess0, LB, UB, optvar,w',consts,rest_fdata);
     
figure
plot(w, rest_fdata,'g');hold on
plot(w, AR1spec_lsq(guess,w,consts),'k')
axis([-4 4 -100 3500])
title('spectrum of the mean')
xlabel('frequency')
legend('data', sprintf('fitted A = %f , rho=%f', guess(1), guess(2)))
guess

function result = AR1spec_lsq(parms , w , consts, data);
A = parms(1);
rho = parms(2);
%B = parms(3);

ar = ones(size(w));
ar = abs(A* ar ./ (1 - rho * exp(-j.*w))) ;
if nargin==4
    result = data-ar;
else
    result = ar;
end
return