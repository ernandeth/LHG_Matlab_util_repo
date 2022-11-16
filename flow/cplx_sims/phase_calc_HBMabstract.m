phis = [];
angs=[];
close all

% Ca = is the arterial volume fraction
for Ca =0.01:0.01: 0.12 ; 

    Gz = 1;  % Flow crusher gradient:  G / cm
    dt = 0.005; % imaging time during which phase is accrued
    vz = 2; % arterial velocity cm/s
    f = 0.015 ; % perfusion rate ml/s/ml
    R1a = 1/1.6;  % arterial R1
    Ttrans = 1.2;  % transit time in seconds
    alpha = 1 * exp(-Ttrans * R1a); %( T1 adjusted inversion efficiency)
    %alpha = 0
    %Ca = 0.10;

    gamma = 2*pi*4.257e3; % rad/s/gauss
    phia = gamma * Gz * vz * dt^2
    phis = [phis phia];

    
    % Cb is the tissue volume fraction
    Cb = 1-Ca;
    % the amount of spins that goes into the tissue is Cf
    Cf = f*(Cb-Ca);
    Cb = 1 - Ca - Cf;
    
    phib = 0;

    % tissue complex signal
    bsig = (Cb - (1-2*alpha)*Cf)*  exp(-i*phib);
    % arterial complex signal
    asig = (1-2*alpha)*Ca*exp(-i*phia);

    % net signal in the control case
    signal1 = ...
        Cb * exp(-i*phib) + ...
        Cf * exp(-i*phib) + ...
        Ca * exp(-i*phia);

    % net signal in the tagged case
    signal2 = ...
        Cb *  exp(-i*phib) + ...
        -Cf *(1-2*alpha) *  exp(-i*phib) + ...
        (1-2*alpha)*Ca*exp(-i*phia);

    % complex subtracted signal:
    asl = signal1 - signal2;
    angs = [angs angle(asl)];
    
    
    % show control and tag vectors
    compass( real(signal1) ,  imag(signal1))
    hold on
    compass( real(signal2) ,  imag(signal2),'r')
    hold off
    drawnow
    pause(0.1)
end

figure
plot( angs)