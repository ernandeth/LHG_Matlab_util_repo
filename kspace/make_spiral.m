function [k_spiral, grads] = make_spiral(dt, shot_duration, fov, phase)
% function [k_spiral, grads] = make_spiral(dt, shot_duration, fov, phase)
% 
% k_spiral  - the trajectory in kspace (complex vector)
% grads - the gradient waveform (complex waveform)
%
% dt        - sample time = 1/BW
% shot_duration     - duration of the spiral
% fov       - field of view
% phase     -  phase of the spiral in radians (use for interleaves)
%


    GAMMA = 26752;  % rad/s/Gauss

    % figure out the fov stuff assuming 1x1cm pixels
    kmax = 1/2; %(samples/ cm)
    
    % spiral gradient waveforms
    t=[0:dt:shot_duration];
    
    %phase = phase*dt/(2*pi);
    
    %t= t + phase;
    t =t.^(2/3);
    spfreq = (pi*fov)/ t(end);
    
    
    k_spiral=  complex(t .* sin(spfreq*t +phase) , t.*cos(spfreq*t+phase));
    k_spiral = k_spiral * kmax/(abs(max(k_spiral))); % normalize
    
    
    grads = diff(k_spiral / dt) * (2*pi)/GAMMA;
    max_grads = max(real(grads));
    

return