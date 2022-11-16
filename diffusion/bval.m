function b=bval(Delta2,delta,Gamp,ramptime)
%function b=bval(Delta2,delta,Gamp,ramptime)
%
%     Delta2 = spacing between lobes (ms)
%
%     delta = width of lobe plateau (ms)
%
%     Gamp = gradient amplitude (Gauss/cm)
%
%     ramptime  (optional) = ramp time for the diffusion gradients (ms).
%     --------------------------------------------------------------------------
%     Assuming ramptime=0 cause a significant difference in the calculated b-value 
%     unless  delta>>ramptime

if nargin<4
    ramptime=0;
end


     gamma=4257.7*2*pi; % (Hz/Gauss)  %NEED THE 2pi to match values reported in the literature

     Gamp=Gamp/10; % (make Gauss/mm)
     Delta2=Delta2/1000; delta=delta/1000; % (ms to s)
     ramptime=ramptime/1000; % (ms to s)

if(ramptime>0)
     delta=delta+ramptime;  
     %%This equation is from Pell et al from MRM 49: 341-350  (2003)%%
     b=(gamma*Gamp*delta).^2*(Delta2-delta/3)+(gamma*Gamp).^2*((ramptime.^3)/30)-(gamma*Gamp).^2*((ramptime.^2)*delta/6);  %in s/mm^2
else
     b=(gamma*Gamp*delta).^2*(Delta2-delta/3);  %in s/mm^2
     b2=2/3*delta*(gamma*Gamp*delta).^2; %Only valid for Delta2=delta. from Haacke eqn 21.34. same result as above for this case
end


bipolar_pair_length_ms=(Delta2+delta+2*ramptime)*1000
GE_forward_spiral_minTE=ceil(bipolar_pair_length_ms)
revspiral_readout_length=26;
GE_reverse_spiral_minTE=ceil(bipolar_pair_length_ms)+revspiral_readout_length

