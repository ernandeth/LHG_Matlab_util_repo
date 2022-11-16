function  obs = gen_signals_vs_200210(parms, delays, aq_parms, dofigs,doSub,b1gr_parms)
% note: defined new OPTIONAL input b1gr_parms containing predefined B1 and
% gr inputs
% specify constants here:
dt = 0.01; % time step for decay/flow calc
dt2 = 0.0001; % time step for rotation calc
rotTime = .005 % time duration of rotation event, assumed to be 5 ms arbitratily
t2 = linspace(0,rotTime,rotTime/dt2);
R2 = 1/0.125;
R2a = 1/0.250;
R1 = 1/1.2;
R1a = 1/1.7;
% Relaxation operator:
R = diag(exp(-dt*[R2, R2, R1, R2a, R2a, R1a] ));
lambda = 0.9;
nslices = 16; % possible this should be in params or just not used? not sure
gam = 4257.746778 * 2 * pi * 1e-2; 
frames = [];
df= 0; % assume off-resonance freq of position 0? not sure if we want this in parms


if (nargin == 0)
    f = 0.01; % 60 ml/min/100g
    bat = 0.12;



    t_tag = 1.5;
    t_delay = 0.2;
    t_adjust = 2;
    t_aq = 0.6 


    labelcontrol = [0 0 1 -1 1 -1];
    Tmax = t_adjust + t_tag + t_delay + t_aq;


    t = linspace(0,Tmax,Tmax/dt);

    B1 = zeros(2,length(t2)); % 2 dimensions by time vector length so that there can be one vector for sel & one for nonsel
    gr = zeros(2,length(t2));

    % perfusion/exchange operator:
    % outflow:
    F = diag(exp(-dt*[f/lambda , f/lambda ,f/lambda , 0,0,0]));
    % inflow:
    F(3,6) = 1-exp(-f*dt/lambda);

    % equilibrium state:
    M0 = [0 0 1 0 0 1]' ;

    % initial state:
    M = [0 0 1 0 0 1]';

    % hi res input function:
    dt2 = 1e-3;   % new temp resoultion
    disp = 0.15; % dispersion
    eff = 0.95;  % inversion efficiency
    tau1 = bat; % leading edge
    tau2 = bat+t_tag; % trailing edge
    aif = vsasl_inputfun(dt2, disp, eff, R1a, tau1, tau2);
    aif = aif(1:dt/dt2:end);
    % use the precomputed input function:
    F(6,6) = 0;
    %%%

    result = ones(Tmax/dt,2); 

    N_frames = 6;
elseif (nargin == 5)
    t_tags = aq_parms.t_tag;  % zero if you have a single pulse
%     t_delays = aq_parms.t_delays;
    t_delays = delays;
    t_adjusts = aq_parms.t_adjusts;
    labelcontrol = aq_parms.labelcontrol;    
    doArtSup = aq_parms.doArtSup;
    ArtSup_delay = aq_parms.ArtSup_delay ; % delay between AS pulse and acqusition
    
    aqwindow = aq_parms.t_aq(1);
    
    f = parms.f;
%     Mtis0 = parms.Mtis0;
    Mtis0 = 1;
    cbva = parms. cbva;
    bat =  parms.bat;
    r1tis =  parms.r1tis;
    flip =  parms.flip;
    alpha_ai = parms.alpha_ai; % arterial inversion efficiency
    alpha_ti = parms.alpha_ti; % tissue inversion efficiency
    alpha_ts = parms.alpha_ts; % T2 weighting in tissue due to arterial suppression 
%     alpha = 0.8;
    Tmax = t_adjust + t_tag + t_delay + t_aq;
    t = linspace(0,Tmax,Tmax/dt);
    
    B1 = zeros(2,length(t)); % 2 dimensions by time vector length so that there can be one vector for sel & one for nonsel
    gr = zeros(2,length(t));
    % perfusion/exchange operator:
    % outflow:
    F = diag(exp(-dt*[f/lambda , f/lambda ,f/lambda , 0,0,0]));
    % inflow:
    F(3,6) = 1-exp(-f*dt/lambda);

    % equilibrium state:
    M0 = [0 0 1 0 0 1]' ;

    % initial state:
    M = [0 0 1 0 0 1]';

    % hi res input function:
    dt2 = 1e-3;   % new temp resoultion
    disp = 0.15; % dispersion
    eff = 0.95;  % inversion efficiency
    tau1 = bat; % leading edge
    tau2 = bat+t_tag; % trailing edge
    aif = vsasl_inputfun(dt2, disp, eff, R1a, tau1, tau2);
    aif = aif(1:dt/dt2:end);
    % use the precomputed input function:
    F(6,6) = 0;
    %%%

    result = ones(Tmax/dt,2); 

    N_frames = length(t_delays);
else
   t_tags = aq_parms.t_tag;  % zero if you have a single pulse
%     t_delays = aq_parms.t_delays;
    t_delays = delays;
    t_adjusts = aq_parms.t_adjusts;
    labelcontrol = aq_parms.labelcontrol;    
    doArtSup = aq_parms.doArtSup;
    ArtSup_delay = aq_parms.ArtSup_delay ; % delay between AS pulse and acqusition
    
    aqwindow = aq_parms.t_aq(1);
    
    f = parms.f;
%     Mtis0 = parms.Mtis0;
    Mtis0 = 1;
    cbva = parms. cbva;
    bat =  parms.bat;
    r1tis =  parms.r1tis;
    flip =  parms.flip;
    alpha_ai = parms.alpha_ai; % arterial inversion efficiency
    alpha_ti = parms.alpha_ti; % tissue inversion efficiency
    alpha_ts = parms.alpha_ts; % T2 weighting in tissue due to arterial suppression 
%     alpha = 0.8;
    Tmax = t_adjust + t_tag + t_delay + t_aq;
    t = linspace(0,Tmax,Tmax/dt);
    
    B1 = b1gr_parms.B1;
    gr = b1gr_parms.gr;   
    % perfusion/exchange operator:
    % outflow:
    F = diag(exp(-dt*[f/lambda , f/lambda ,f/lambda , 0,0,0]));
    % inflow:
    F(3,6) = 1-exp(-f*dt/lambda);

    % equilibrium state:
    M0 = [0 0 1 0 0 1]' ;

    % initial state:
    M = [0 0 1 0 0 1]';

    % hi res input function:
    dt2 = 1e-3;   % new temp resoultion
    disp = 0.15; % dispersion
    eff = 0.95;  % inversion efficiency
    tau1 = bat; % leading edge
    tau2 = bat+t_tag; % trailing edge
    aif = vsasl_inputfun(dt2, disp, eff, R1a, tau1, tau2);
    aif = aif(1:dt/dt2:end);
    % use the precomputed input function:
    F(6,6) = 0;
    %%%

    result = ones(Tmax/dt,2); 

    N_frames = length(t_delays);
end

for i = 1:N_frames
    for n = 1:(t_adjust/dt)
        M = (F * R) * M + (eye(6) - F*R) * M0;
        result(n,1) = M(3,1);
        result(n,2) = M(6,1);
    end

    m=1;
    M(3) = 1-2*eff;  % invert the stationary spins.
    % note: to be replaced by B1/gr-caused rotation 
%   for j = 1:length(t2)
%         b1x = real(B1(1,j));
%         b1y = imag(B1(1,j));
%         Bz = 2 * pi * df / (gam * 100) % + gr(1,j) * position?  need to define spatial position wrt gradient
%         % incorporate off-resonance frequency into Bz?
%         B_abs = sqrt(b1x^2 + b1y^2 + Bz^2);
%         theta = gam * 1e2 * B_abs * dt;
%         if B_abs > 0
%                 ux = b1x / B_abs;
%                 uy = b1y / B_abs;
%                 uz = Bz  / B_abs;
%             else
%                 ux = 0;
%                 uy = 0;
%                 uz = 0;
%             end
% 
%             % Rotations are left-handed
%             c = cos(-theta);
%             s = sin(-theta);
%             one_minus_c = 1 - c;
% 
%             R11  = c + ux * ux * one_minus_c;
%             R22  = c + uy * uy * one_minus_c;
%             R33  = c + uz * uz * one_minus_c;
%             R12a = ux * uy * one_minus_c;
%             R12b = uz * s;
%             R13a = ux * uz * one_minus_c;
%             R13b = uy * s;
%             R23a = uy * uz * one_minus_c;
%             R23b = ux * s;
% 
%             Rot = [R11        , R12a - R12b, R13a + R13b, 0,0,0;
%                  R12a + R12b, R22        , R23a - R23b, 0,0,0;
%                  R13a - R13b, R23a + R23b, R33        , 0,0,0;
%                  0          , 0          , 0          , R11        , R12a - R12b, R13a + R13b
%                  0          , 0          , 0          ,R12a + R12b, R22        , R23a - R23b;
%                  0          , 0          , 0          ,R13a - R13b, R23a + R23b, R33];
%              M = Rot * M;
%       end
    for n = (t_adjust/dt):(t_adjust+t_tag+t_delay)/dt-1
        M = (F * R) * M + (eye(6) - F*R) * M0; 

        if labelcontrol(i)==0
            if m<length(aif)
                M(6) = aif(m); m=m+1;
            end
        end
        result(n,1) = M(3,1);
        result(n,2) = M(6,1);
    end

    for n = (t_adjust + t_tag+t_delay)/dt : Tmax/dt
        M = (F * R) * M + (eye(6) - F*R) * M0; 
        M(3) = M(3) * (3)*cos(flip)^nslices; % note: also replace with rotation loop following this section?
        result(n,1) = M(3);
        result(n,2) = M(6);
    end
    frames(i) = M(3);
end


return