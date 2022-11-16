%% Calculate b-value for pair of symmetric bi-polar trapezoid gradients
%
% b = bValueA(Gmax, f, r, sep, gamrad)
%
% in:
%      Gmax   - maximum gradient amplitude (units/cm)
%      f      - flat top duration (ms)
%      r      - ramp time (ms)
%      sep    - seperation between the centre of the trapezoids (ms)
%      gamrad - gyromagnetic ratio (rad/s/units)
%
% out:
%      b - calculated diffusion b-value (s/mm^2)
%
% example usage:
%
% b = bValueA(1 G/cm, 2 ms, 0.5 ms, 10 ms, gyroratio('Hz/G')*2*pi)
%
% Reference: Bernstein et al., Handbook of MRI Pulse Sequences (2004) page 278
%
% Written by Joseph G. Woods, CFMRI, UCSD, April 2020

function b = bvalueA(Gmax, f, r, sep, gamrad)

Gmax_mm = Gmax / 10;  % convert from units/cm to units/mm
f_s     = f   * 1e-3; % convert from ms to s
r_s     = r   * 1e-3; % convert from ms to s
sep_s   = sep * 1e-3; % convert from ms to s

delta = f_s + r_s;

b = gamrad^2 * Gmax_mm^2 * ( delta^2*(sep_s-delta/3) + r_s^3/30 - delta*r_s^2/6 );

end
