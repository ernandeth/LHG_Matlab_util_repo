function [inputChirp, tsamp] = getInputChirp(procpar)
% function [inputChirp] = getInputChirp(procpar)
%
%  Script to get input chirp from files
%  the parameters defining the specific chirp were
%  stored in the procpar file

% Get constants
chirpAmp = -0.5;
f0 = procpar.f0;
f1 = procpar.f1;
at = procpar.at;
dt = 4e-6;          % Resolution
t = 0: dt :(at-dt);

SR = procpar.gmax/procpar.trise;
dt = 1/procpar.sw;
tsamp = 0: dt : at-dt;
gamma = 4257;  

%t = tsamp;

for uu = 1:numel(t)
    
    f = f0 + (f1-f0)*t(uu)/at/2;
    y = sin(2*pi*f*t(uu));
    
    % Figure out amplitude
    var1 = 2*pi*(f0*(1-t(uu)/at)+f1*t(uu)/at);
    
    
     if (SR/var1<= chirpAmp)
          Amp = SR/var1;
     else
         Amp = chirpAmp;
     end
    inputChirp(uu) = Amp*y;
    
end

% resample input chirp 
inputChirp = interp1(t,inputChirp,tsamp, 'sinc');


end
