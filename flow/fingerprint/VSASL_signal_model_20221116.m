inp_fun = ones(Nsamp,1);
prep1_fun = zeros(Nsamp,1);
prep2_fun = zeros(Nsamp,1);
readout_fun = zeros(Nsamp,1);


% possible events:

% Prep pulse 1:  labeling pulses 
%----------------------------------
% effect on tissue: (inversion efficiency)
a1t_sel = 0.8;
a1t_nonsel = 0.2;
% effect on blood: (inversion efficiency)
a1a_sel = 0.9
a1a_nonsel = 0.1;


% Prep pulse 2:  arterial suppression pulses 
%------------------------------------------
% effect on tissue: (inversion efficiency)
a2t_sel = 0.1;
a2t_nonsel = 0.1;
% effect on blood: (inversion efficiency)
a2a_sel = 0.5
a2a_nonsel = 0.1;

% FSE readout
%---------------------------------------
% effects of readout on the magnetization
a3t = 0.5;
a3a = 0.5;


switch prep1_fun(n)
    case 1
        
    case 0
        
end

switch prep2_fun(n)
    case 1
        
    case 0
end
