phis = [];
angs=[];
close all

% Arterial Volume of a ROI normalized to 100g ~ 15% of 4% total CBV
% Arteriolar volume of a ROI normalized to 100g ~ 10% of 4% total CBV
% Capillary volume of a ROI normalized to 100g ~ 10% of 4% total CBV
% Venule volume of a ROI normalized to 100g ~ 20% of 4% total CBV
% Venous volume of a ROI normalized to 100g ~ 45% of 4% total CBV
Ca = 0.20;

% Ca = is the arterial volume fraction
f = 0.015;

% Cb is the tissue volume fraction
Cb = 1-Ca;

Gz = 4;  % Flow crusher gradient:  G / cm
dt = 0.002; % imaging time during which phase is accrued
Delta = 0.0088;



alpha = 0.9
R1a = 1/1.6;  % arterial R1
Ttrans = 1.2;  % transit time in seconds
alpha = alpha*exp(-Ttrans * R1a); %( T1 adjusted inversion efficiency)

%alpha = 0

% reduce the arterial contribution gradually now
for vz =5:-0.2: 0.1 ;
    % fraction of the signal in arteries (not capillaries or tissue)
    %Ca = 0.5;

    %vz = 2;
    % arterial velocity cm/s
    % f = 0.015 ; % perfusion rate ml/s/ml
    % using the relationship from our model
    % vz = 2 + (0.2)*(f-0.01)*(2/0.01)

    gamma = 2*pi*4.257e3; % rad/s/gauss
    % phia = gamma * Gz * vz * (-2*dt)^2;
    phia = -gamma * Gz * vz * Delta* dt; % see notebook p. 53 (2.21.08)
    phis = [phis phia];


    % the amount of spins that goes into the tissue is Cf
    [Ca Cb ];
    phib = 0;

    % tissue complex signal
    bsig = Cb *  exp(-i*phib);
    % arterial complex signal
    asig = (1-2*alpha)*Ca*exp(-i*phia);

    % net signal in the control case
    signal1 = ...
        Cb * exp(-i*phib) ...
        + Ca * exp(-i*phia);

    % net signal in the tagged case
    signal2 = ...
        Cb *  exp(-i*phib) ...
        + (1-2*alpha)*Ca*exp(-i*phia);

    % complex subtracted signal:
    asl = signal1 - signal2
    
    [phia phib angle(signal1) angle(signal2) angle(asl)]
    angs = [angs angle(signal1)-angle(signal2)];
    
    
    
    % show controla and tag 
    compass( real(signal1) ,  imag(signal1))
    hold on
    compass( real(signal2) ,  imag(signal2),'r')
    compass( real(asl), imag(asl),'g')
    hold off
    drawnow
    pause(0.1)
end
VENC = 2*pi / (gamma * Gz * Delta* dt) % see notebook p. 53 (2.21.08)
figure
plot( angs)

% some calculations from  Grubb
% flow increase of 20%
% (v/v0) = (f/f0)^0.38

relV = (1.3)^0.38

