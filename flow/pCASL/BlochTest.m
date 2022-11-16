close all
dt = 1e-3; % 1 us expressed in millliseconds
nstep = 50 / dt;  % 50 ms @  1us steps.

T1 = 1200;
T2 = 50;

Mi = [0, 0, 1]; % initial condition: tipped into x axis.

B0 = 3;    % keep things in Tesla

% gradients in Gauss/cm.
Gx = zeros(nstep,1);
Gy = zeros(nstep,1);
Gz = zeros(nstep,1);

Bx = zeros(nstep,1);
By = zeros(nstep,1);
Bz = zeros(nstep,1);

Gx(20/dt:30/dt) = 1;
Gx(31/dt:40/dt) = -1;

% position in cm
dx = 1;
dy = 0;
dz = 0;

% Bz field caused by gradients expressed in Tesla
Bz = (dx*Gx + dy*Gy + dz*Gz) * 1e-4;

% B1 field in mGauss
%define a 6 milllisecond sinc pulse with two lobes on each side
rfdt = 4*pi / (6/dt) ;
mySinc = sinc(-2*pi : rfdt : 2*pi);
% scale the pulse to 120mGauss max
mySinc = mySinc * 120;  
% conversion from mGauss to Tesla
mySinc = mySinc * 1e-7;

RFpulse = mySinc;

% shift the frequency of the RF pulse by modulating the RF
% add phase to the pulse by adding phase to the carrier
RF_dw =0 ;
RF_theta = pi/8;

FMx = cos(RF_dw  * 2*pi * [1:length(mySinc)] * dt - RF_theta);
FMy = sin(RF_dw  * 2*pi * [1:length(mySinc)] * dt - RF_theta);

myRF_x = mySinc .* FMx;
myRF_y = mySinc .* FMy;

% integrate the pulses into the sequence
Bx(1100:1100+length(RFpulse)-1) = myRF_x;
By(1100:1100+length(RFpulse)-1) = myRF_y;

beff = [Bx By  Bz];

M = blochsim(Mi, beff, T1, T2, dt, nstep);
figure
subplot(411) , plot(Gx)
subplot(412) , plot(Gy)
subplot(413) , plot(Gz)
subplot(414) , plot(Bx), hold on, plot(By,'r')

figure
subplot(311) , plot(M(:,1))
subplot(312) , plot(M(:,2))
subplot(313) , plot(M(:,3))
