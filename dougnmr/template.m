% all times in ms, all freqs in kHz, all distances in cm, all B's in T

% load object; % loads in file of objects

gambar = 42570;               % gamma/2pi in kHz/T
gam = gambar*2*pi;            % gamma/2pi in kHz/T

% simulation values
dt = 16e-3;                   % ms
te = 12.0;                    % ms
endtime = 17;                 % ms
time = [0:dt:endtime]';       % ms
npts = length(time);          % number of time points for simulation
m = zeros([npts 3 ]);         % place holders
beff = zeros([npts 3 ]);      % place holders

% define 90 RF
pwrf90 = 1.6;                 % in ms
sincper = pwrf90/4;           % in ms
rfsteps = pwrf90/dt;
rftime = [-(rfsteps-1)/2:(rfsteps-1)/2]'.*dt;
rfshape = hanning(rfsteps).*sinc(rftime./sincper);
a_rf90 = ???                  % in T
b1_90 = a_rf90.*[rfshape; zeros([npts-rfsteps 1])]; % in T

% define gz
% slthick = 1;                  % in cm
% a_gz1 = ??? % in T/cm
% pwgz1 = pwrf90;
% a_gz2 = -a_gz1;
% pwgz2 = pwrf90/2;
% gz =  (time < pwgz1) .* a_gz1 ...
%        + (time >= pwgz1).*(time < (pwgz1+pwgz2)) .* a_gz2;

% define 180 RF
pwrf180 = 0.8;                 % in ms
n180 = round(pwrf180/dt);
a_rf180 = ???                  % in T
nloc180 = round((te/2 + pwrf90/2 - pwrf180/2)/dt);
b1_180 = a_rf180.*[ zeros(nloc180,1); ones(n180,1); ...
		    zeros(npts-nloc180-n180,1)]; % in T

% define gx

% define gy

% for slice = 1:2  % slice loop
    % 
    % modify B1 here for other slices
    % 

    % for pe = 1:npe % phase encode loop
        %
	% assign a_gy values here
	%
        
        % for obj = 1:nobj % object loop
	    %
	    % set beff here in terms of RF, gradients, location of spin
	    % beff should be of size npts x 3 and have units of Tesla
	    beff = ???
	    %
	    % also set up T1, T2 for each object
	    T1 = ??? 
            T2 = ???
	   
	    m0 = [0 0 1];  % initial magnetization
            m = blochsim(m0,beff,T1,T2,dt,npts);

	    % defined received (sampled) signal from m
            % xpos = [-nread/2+1:nread/2]/nread*30;
            % plot(xpos,abs(ift(sig)));
        % end  % object loop

	%
	% place each phase encode line into the acquisition matrix
	% specifically, load M(kx,ky)
	%

    % end  % phase encode loop
    %
    % reconstruct the image here
    % show(abs(image));
    % label and print
    %

% end  % slice loop

