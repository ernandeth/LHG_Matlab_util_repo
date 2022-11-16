function pulse_out = scaleB1(pulse_in, flip);
%
% function pulseout = scaleB1(pulse_in, flip);
% pulse_in is the B1 pulse waveform.  Each point is 1 usec.
% flip is in degrees
%
% pulse_output is the scaled version in Gauss
% 
% example:  the GE scanner calibration pulse
% s1 = sinc(linspace(-2,2,3200)) .* hanning(3200)'; plot(s1)
% s1_scaled = scaleB1(s1,180); max(s1_scaled)
%

GAMMA = 26752;  % rad/s/Gauss
pulse_in = pulse_in/max(pulse_in);
dt = 1e-6;
B1area = sum(pulse_in)*dt;
pulse_out = pulse_in * deg2rad(flip) / (B1area * GAMMA);
return



s1 = sinc(linspace(-4,4,6400)) ; plot(s1)

s1 = sinc(linspace(-8,8,6400)) .* hamming(6400)'; plot(s1)
s1_scaled = scaleB1(s1,180); max(s1_scaled)
