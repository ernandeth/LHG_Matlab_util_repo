function result = MTspec_lsq( guess, dw,  parms, data)
%function MTspec_lsq( [ka_guess f_guess], dw, parms {, Data})

gamma = 4.258;  % mGauss / Hz

H1 = parms(1);  % mGauss
T1 = parms(2);	% sec. 
T2 = parms(3);	% sec
f = parms(4);   % no units
T1b = parms(5); % sec
T2b = parms(6); %sec

w1 = gamma * H1;
R1 = 1/T1;
R2 = 1/T2;
R1b = 1/T1b;
R2b = 1/T2b;

ka = guess(1);  % sec^-1
%f = guess(2);

mz = -R1./(R1 - ka*(1-f) - ((gamma*H1)^2 * T2)./(1 + (dw*T2)));

% Henkelmann's equation (MRM 29,p.759, 1993):
% Mz is expressed relative to Mz(0)
Rrfa = w1^2 * T2  ./ (1 + (2*pi*dw * T2).^2);
Rrfb = w1^2 * T2b ./ (1 + (2*pi*dw * T2b).^2);
%keyboard 
mz = (R1b*ka*f + Rrfb*R1 + R1b*R1 + R1*ka )./ ...
   ( (R1 + Rrfa + ka*f) .* (R1b + Rrfb + ka) - ka*ka*f);

if nargin>3
    result = mz-data;
else
    result = mz;
end

return
