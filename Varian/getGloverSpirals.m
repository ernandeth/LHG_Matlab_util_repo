%  Script to create and export spirals using Gary Glover's Simple Analytic
%  Spirals spiral waveform (1999, Mangeit Rsonance in Medicine).  In it's finished form, the script will export
%  Gx and Gy as well as Kx and Ky with the only inputs being a procpar file
%  containing all the necessary parameters.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% S P E C I A L  N O T E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The 205 gradient coils have a 44 and 48 microsecond delay between when a gradient is started and when it is actually played
% For the X and Y Gradient directions respectively.
% This delay is hard coded into this .m file so that the spirals used for regridding and reconstruction
% are as identical as possible to the spirals actually played out.
%
% We have also noticed that output gradients are on averate about 95% of
% the input waveform.  So that is entered in here too.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Gx Gy] = getGloverSpirals(procpar);
% Set up contants
q = 5;   % Number must be 1 or larger.  See referenced paper
gamma = 4257.707747;      % Gyromagnetic ratio
S0 = procpar.spSlewRate;   % Max Slewrate of the gradients
fov = procpar.lro;        % Field of view cm
Nint = procpar.Nint;     %# of interleaves (also called shots)
Amp = Nint/(2*pi*fov);   %Amplitude of spiral k-space trajectory
beta = S0*gamma/Amp;    % Gradient Amplification factor
a2 = (9*beta/4)^(0.333334);  % Useful term to compute
spGmax = procpar.Gamp;    % Max gradient in G/cm
mat = procpar.Mat;         % desired reconstruction matrix size.

% Useful estimates of the timeit takes to acquire the spiral

Ts = 0.666667/Nint*sqrt(((pi*mat)^3)/(2*pi*gamma*fov*S0)); % Time to completion of the slew rate limited case
ts = (3*gamma*spGmax/(2*Amp*a2*a2))^3;  % Time to maximum gardient amplitude acheived.

% Slew Rate limited Spiral calculation
i = 1;  % keep track of the numberof points in the spiral
sp_thetaX = 0;    % Angle orotTrajf the k-space trajectory
sp_thetaY = 0;    % Angle of the k-space trajectory
% t = 0;
tX = -0.000045;          % keep triack of time, initialize at t = -44 us.  See special note above;
tY = -0.000048;            % keep triack of time, initialize at t = -44 us.  See special note above;
while ( sp_thetaX <= mat*pi*Nint)
    if tX<0 && tY <0
        Gx(i) = 0;
        Gy(i) = 0;
    elseif tY <0 && tX >=-1e-9
        xX = tX^1.333333;
        yX = q+0.5*beta/a2*xX;
        sp_thetaX = 0.5*beta*tX*tX/yX;
        dthdtX = beta*tX*(q+0.166667*beta/a2*xX)/(yX*yX);
        cX = cos(sp_thetaX);
        sX = sin(sp_thetaX);
        Gx(i) = 1*Nint/(fov*gamma*2*pi)*dthdtX*(cX-sp_thetaX*sX);
        Gy(i) = 0;
    else
        xX = tX^1.333333;
        yX = q+0.5*beta/a2*xX;
        sp_thetaX = 0.5*beta*tX*tX/yX;
        dthdtX = beta*tX*(q+0.166667*beta/a2*xX)/(yX*yX);
        cX = cos(sp_thetaX);
        sX = sin(sp_thetaX);
        Gx(i) = 1*Nint/(fov*gamma*2*pi)*dthdtX*(cX-sp_thetaX*sX);
        xY = tY^1.333333;
        yY = q+0.5*beta/a2*xY;
        sp_thetaY = 0.5*beta*tY*tY/yY;
        dthdtY = beta*tY*(q+0.166667*beta/a2*xY)/(yY*yY);
        cY = cos(sp_thetaY);
        sY = sin(sp_thetaY);
        Gy(i) = 1*Nint/(fov*gamma*2*pi)*dthdtY*(sY+sp_thetaY*cY);
        
        % if statement to ensure taht maximum gradient amplitude is not exceeded
        Gabs = hypot(Gx(i), Gy(i));
        if (Gabs>= spGmax)
            if (Gabs > hypot(spGmax/sp_thetaX,spGmax))
                i = i+1;
                break;
            end
        end
    end
    i = i+1;
    tX = tX+procpar.spResolution;
    tY = tY+procpar.spResolution;
end

% Check to see if spiral sequence in finished.  If not, use amplituded limited design
thetas = sp_thetaX; ts = tX;
if thetas < (mat*pi/Nint)
    Tta = ts+pi*Amp/(2*pi*gamma)/spGmax*(mat*pi/Nint)^2 - (thetas*thetas);
    t = ts+procpar.spResolution ; % Initialize t again.
    
    % Compute the amplitude limited waveform
    sp_theta = sp_thetaX;
    while sp_theta <= (mat*pi/Nint)
        sp_theta = sqrt(thetas*thetas+2*gamma*spGmax/Amp*(t-ts));
        dthdt = gamma*spGmax/Amp/sp_theta;
        c = cos(sp_theta);
        s = sin(sp_theta);
        Gx(i) = 1*spGmax*(c/sp_theta-s);
        Gy(i) = 1*spGmax*(s/sp_theta+c);
        i = i+1;
        t = t+procpar.spResolution;
    end
else
    t = (tX+tY)/2;
end

% Ramp down the gradients to zero in 75 points at 4 us each

j = i-1;  % Initialize ramping parameters
n = 0;   % initialize ramping parameters
while n <=75
    Gx(i) = Gx(j)*(1-n/75.0);
    Gy(i) = Gy(j)*(1-n/75.0);
    i = i+1;
    n = n+1;
end

% Add one more zero point to match the spirals computed in spiralsgl.c
Gx = Gx';
Gy = Gy';
t = t+procpar.spResolution;
end

