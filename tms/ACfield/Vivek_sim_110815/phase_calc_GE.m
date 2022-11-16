function [phase_vec,sin_wave_vec] = phase_calc_GE(procpar)
%function [phase vector, sine wave] = phase_calc_GE(procpar)
%
%Accepts scan information data 'procpar' as input . Procpar should have the
%members size of image,acqusition time, freq of sinusoidal phase and starting phase .
%Generates the sine wave and calculates its integral using the Simpson's
%rule
%
%VRB U of M 05/24/2010 | bhatiavr@umich.edu


%----------------------Phase calcs for Bfield------------------------------
%  sin_func = @(t) exp(w*1i*t + intial_phase);
% time = linspace(0,(procpar.at*procpar.nv), (procpar.np/2)^2);


time = linspace(0,(procpar.at),(procpar.np/2)); %Sampling time points 
phase_vec = zeros(1,procpar.np/2);
dt = (time(2)-time(1)); %Sampling time points
freq = procpar.freq;
initial_phase = procpar.initial_phase; % in degrees


% %--------------------------------------------------------------------------
% %Sinusoidal wave for each data acquistion phase 
% sin_wave = sinewave(time,freq,initial_phase,1);
% 
% 
% for itime = 1:procpar.np/2
%     phase_vec(1,itime) = quad(@(t)sinewave(t,freq,initial_phase,2),time(1),time(itime)); % Calclate integral of sinewave by Simpsons method
% end
% 
% % phase_vec(1,:) = cumsum(imag(sin_wave)) .* dt;
% % sin_wave = repmat(sin_wave,1,procpar.nv);
% phase_vec = repmat(phase_vec,1,procpar.nv);
% 
% 
% % eulers = @(t,freq,initial_phase) exp(1i .* (2.*pi.*freq.*t + initial_phase));
% % imag_eulers = @(t) imag(exp(1i .* (2.*pi.*freq.*t + initial_phase)));
% 
% 
% % 
% % count = 2;
% % for ii  = 1:procpar.nv
% % %     subplot(211)
% %     plot(phase_vec((1:1:procpar.np/2)+(procpar.np/2)*(count-2),1));
% % %     subplot(212)
% % %     plot(time,imag(sin_wave(1,(1:1:procpar.np/2)+(procpar.np/2)*(count-2))),'g');
% % 
% %     %     plot(time,cumsum(imag(sin_wave(1,(1:1:procpar.np/2)+(procpar.np/2)*(count-2))))/32, 'g');
% %     pause;
% %     clf;
% %     count = count + 1;
% % end
% % 
% 
% sin_wave_vec = sin_wave;

%--------------------------------------------------------------------------


%Sinusoidal wave for each data acquistion phase 
sinewave=@(t) sin(2*pi*freq*t + initial_phase*pi/180);


for itime = 1:procpar.np/2
    phase_vec(1,itime) = quad(@(time)sinewave(time),time(1),time(itime)); % Calclate integral of sinewave by Simpsons method
end

% phase_vec(1,:) = cumsum(sinewave(time)) .* dt; %This was not giving a good
                                                 % fit
% sin_wave_vec = repmat(sinewave(time),1,procpar.nv);

phase_vec = repmat(phase_vec,1,procpar.nv); % Sine wave is repeated for each phase encode



% % Display sine wave and its integral for each data acqustion phase
% count = 2;
% for ii  = 1:procpar.nv
% %     subplot(211)
%     plot(phase_vec(1,(1:1:procpar.np/2)+(procpar.np/2)*(count-2)));
% %     subplot(212)
% %     plot(time,sinewave(time),'g');
%     pause;
%     clf;
%     count = count + 1;
% end
% 
sin_wave_vec = sinewave(time);


return
