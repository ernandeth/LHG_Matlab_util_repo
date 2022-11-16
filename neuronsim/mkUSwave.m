function wv = mkUSwave(cf, duty, prf, stimlength, duration)
% function wv = mkUSwave(cf, duty, prf, stimlength, duration)
% the time resolution is 0.1 microsecond.
%
%{
cf = 0.5*1e6;   % carrier frequency in Hz
duty = 0.3;     % duty cycle of PRF wave
prf = 1500;     % pulse repetition frequency in Hz
stimlength= 200e-3; % duration of the stimulus within the whole time interval in seconds
duration = 1;   % duration of the time interval in seconds;
%}

dt = 1e-7;  % time units are 0.1 microseconds

t = linspace(0, duration, duration/dt);
wv_stim = zeros(size(t));
wv_prf = zeros(size(t));

wv_carrier = sin(2*pi*cf*t);

wv_stim(1:stimlength/dt) = 1;

cyc_length = 1/prf/dt;  % duration of the whole burst cycle
tbd = duty/prf/dt;      % Tone Burst Duration

for n=1:cyc_length:length(wv_prf)
    wv_prf(n:n+tbd) = 1;
end

wv = wv_carrier .* wv_stim .* wv_prf;

plot(t,wv)

return
