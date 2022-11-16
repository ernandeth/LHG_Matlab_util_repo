tau = 1
iti = 16

corr_ana = bold_shift_ana10(0.3, tau, iti);

save tau=1.mat corr_ana
xs = 0:100:2000;
ys = 100:12:460;

%tc = [3:2:10, 15:10:30, 40:100:300];
%tc = [1:5:100 ];

%%%%%%%%%

clear
tau = 3
iti = 16

corr_ana =  bold_shift_ana10(0.3, tau, iti);

save tau=3.mat corr_ana

%%%%%%%%%

clear
tau = 5
iti = 16

corr_ana =  bold_shift_ana10(0.3, tau, iti);

save tau=5.mat corr_ana

%%%%%%%%%

clear
tau = 6
iti = 16

corr_ana =  bold_shift_ana10(0.3, tau, iti);

save tau=6.mat corr_ana

%%%%%%%%%

clear
tau = 8
iti = 16

corr_ana =  bold_shift_ana10(0.3, tau, iti);

save tau=8.mat corr_ana
%%%%%%%%%

clear
tau = 10
iti = 16

corr_ana =  bold_shift_ana10(0.3, tau, iti);

save tau=10.mat corr_ana

%%%%%%%%%

clear
tau = 6
iti = 2

corr_ana =  bold_shift_ana10(0.3, tau, iti);

save iti=02.mat corr_ana

%%%%%%%%%

clear
tau = 6
iti = 5

corr_ana =  bold_shift_ana10(0.3, tau, iti);
save iti=05.mat corr_ana

%%%%%%%%%

clear
tau = 6
iti = 10

corr_ana =  bold_shift_ana10(0.3, tau, iti);

save iti=10.mat corr_ana

%%%%%%%%%

clear
tau = 6
iti = 15

corr_ana =  bold_shift_ana10(0.3, tau, iti);

save iti=15.mat corr_ana
%%%%%%%%%

clear
tau = 6
iti = 30

corr_ana =  bold_shift_ana10(0.3, tau, iti);

save iti=30.mat corr_ana

%%%%%%%%%

clear
tau = 6
iti = 45

corr_ana =  bold_shift_ana10(0.3, tau, iti);

save iti=45.mat corr_ana
%%%%%%%%%

%  ....Plots ....


xs = -1: 0.1:1;
ys = 0.3:0.3:3;

%tc = [3:2:10, 15:10:30, 40:100:300];
tc = [1:5:100 ];


load iti=02
%subplot 231, [c h] = contour(xs,ys ,corr_ana), title('Mean ITI = 2 sec'),...
%   xlabel('Time Shift'), ylabel('SNR') ;
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 231, plot(xs,  t,  'k'), title('Mean ITI = 2 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8,'fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25]) ;


load iti=05
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 232, plot(xs,  t, 'k'), title('Mean ITI = 5 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25]) ;
%subplot 232, [c h] = contour(xs,ys ,corr_ana), title('Mean ITI = 5 sec'),...
%   xlabel('Time Shift'), label('t score'), axis([-1 1 0 25]), clabel(c,h,[1:5:26]);




load iti=10
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 233,  plot(xs, t, 'k'), title('Mean ITI = 10 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25])  ;
%subplot 233,[c h] = contour(xs,ys ,corr_ana), title('Mean ITI = 10 sec'),...
%   xlabel('Time Shift'), ylabel('SNR'), clabel(c,h,[1:5:26]);



load iti=15
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 234,  plot(xs, t, 'k'), title('Mean ITI = 15 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25]);
%subplot 234,[c h] = contour(xs,ys ,corr_ana), title('Mean ITI = 15 sec'),...
%   xlabel('Time Shift'), ylabel('SNR'), clabel(c,h,[1:5:26]);



load iti=30
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 235,  plot(xs, t, 'k'), title('Mean ITI = 30 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25]) ;
%subplot 235,[c h] = contour(xs,ys ,corr_ana), title('Mean ITI = 30 sec'),...
%   xlabel('Time Shift'), ylabel('SNR'), clabel(c,h,[1:5:26]);


load iti=45
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 236,  plot(xs, t, 'k'), title('Mean ITI = 45 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25]) ;

%subplot 236,[c h] = contour(xs,ys ,corr_ana), title('Mean ITI = 45 sec'),...
%   xlabel('Time Shift'), ylabel('SNR'), clabel(c,h,[1:5:26]);

% ....More Plots ....
figure

load tau=1
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 231,  plot(xs, t, 'k'), title('Tau = 1 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25]) ;
%subplot 231,[c h] = contour(xs,ys ,corr_ana), title('Tau = 0.10 sec'),...
%   xlabel('Time Shift'), ylabel('SNR'), clabel(c,h,[1:5:26]);

load tau=3
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 232,  plot(xs, t, 'k'), title('Tau = 3 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25]) ;
%subplot 232,[c h] = contour(xs,ys ,corr_ana), title('Tau = 0.50 sec'),...
%   xlabel('Time Shift'), ylabel('SNR'), clabel(c,h,[1:5:26]);

load tau=5
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 233,  plot(xs, t, 'k'), title('Tau = 5 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25]) ;
%subplot 233,[c h] = contour(xs,ys ,corr_ana), title('Tau = 1.0 sec'),...
%   xlabel('Time Shift'), ylabel('SNR'), clabel(c,h,[1:5:26]);

load tau=6
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 234,  plot(xs, t, 'k'), title('Tau = 6 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25]) ;
%subplot 234,[c h] = contour(xs,ys ,corr_ana), title('Tau = 1.5 sec'),...
%   xlabel('Time Shift'), ylabel('SNR'), clabel(c,h,[1:5:26]);

load tau=8
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 235,  plot(xs, t, 'k'), title('Tau = 8 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25]) ;
%subplot 235,[c h] = contour(xs,ys ,corr_ana), title('Tau = 2.0 sec'),...
%   xlabel('Time Shift'), ylabel('SNR'), clabel(c,h,[1:5:26]);

load tau=10
t=corr_ana; %*sqrt(1198) ./ sqrt(1 - corr_ana.^2);
subplot 236,  plot(xs, t, 'k'), title('Tau = 10 sec','fontsize',8),...
   xlabel('Time Shift','fontsize',8), ylabel('t score','fontsize',8), axis([-1 1 0 25]) ;
%subplot 236,[c h] = contour(xs,ys ,corr_ana), title('Tau = 2.5 sec'),...
%   xlabel('Time Shift'), ylabel('SNR'), clabel(c,h,[1:5:26]);

