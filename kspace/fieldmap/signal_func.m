% function result = signal_func(estimate,t,n,m,xx,yy,kx,ky,data)
function result = signal_func(estimate,parms,data1)

% disp('func called');

global diagn count_iter % Define global variables


% Assigning important parameters from input data
m = parms.m;
n = parms.size;
xx = parms.xcord;
yy=parms.ycord;
kx = parms.kx;
ky = parms.ky;
time = parms.t;
gambar = parms.gambar;


% Memory allocation
signal_estimate = zeros(1,n*n);
Bfield_estimate_AC = zeros(n,n,n*n);


Bfield_amp = estimate; % Estimate
freq = 10e3; % We know frequency
phase = 0; % We know phase

signal = data1; % Calculated signal


bf_func = @(t) Bfield_amp .* cos(2*pi*freq*t + phase) ; % Model for AC magnetic field


Bfield_estimate_AC(:,:,1) = bf_func(time(1));
Bf_AC = Bfield_estimate_AC(:,:,1);
Bf_AC = Bf_AC(:);

for r = 1 :  n*n
    
    signal_estimate(1) = (m(r) * exp(-1i * 2 * pi * (kx(1) * xx(r) + ky(1) * yy(r))) * exp(-1i * gambar * Bf_AC(r) ))+ signal_estimate(1);
    
end


for itime = 2: n*n
    
    %-----------------------------------Calculation for Bfield-------------
    
    Bfield_estimate_AC(:,:,itime) = (bf_func(time(itime))) + Bfield_estimate_AC(:,:,itime-1);
    Bf_AC = Bfield_estimate_AC(:,:,itime);
    Bf_AC = Bf_AC(:);
    
    
    %----------------------------------------------------------------------
    
    for r = 1 :  n*n
        signal_estimate(itime) = (m(r) * exp(-1i * 2 * pi *((kx(itime)*xx(r)) + (ky(itime)*yy(r)))) * exp(-1i * gambar * Bf_AC(r))) + signal_estimate(itime);
    end
    
end



result =abs(signal -  signal_estimate); % Objective function


count_iter=count_iter+1;

% Save important results
diagn(count_iter).Bf = estimate;
diagn(count_iter).res = sum(result);
diagn(count_iter).sig = signal_estimate;


end



