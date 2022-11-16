% simulate field
vdim = 1000; % microns
vvol = vdim^3;
mdim = 64;

% make a field to put stuff in.
ar = [-mdim/2:mdim/2-1]./mdim.*vdim;
[xx,yy,zz]=ndgrid(ar,ar,ar);i

% define constants 
B0 = 3;
gambar = 42.58e6;
f0 = B0*gambar;  % resonant frequency
chi = -9.6e-6;   % Hz

% arrays of bubble radii and number of bubbles
radii =    [100 50 25 10 100  50 25 ]; 
nBubbles = [ 8   8  8  8  1   8  64 ]; 

simn=3
%for simn = 1:length(radii)
	nbub = nBubbles(simn);
	radius = radii(simn);
	bvol = 4/3*pi*radius^3;
	vfrac = nbub*bvol/vvol;

	for trial = 1:10
		% bubble position
		pos = (rand([nbub 3])-.5)*vdim;

		% define a delta Frequency map
		% C reflects the ratio between the susceptibilites
		% in and outside the sphere
		df = zeros(size(xx));
		C = chi./(3+chi).*f0;

		% loop through the bubbles (theis effect is superimposed)
		for lp = 1:nbub
			% make a field of distances to the bubble position
			rr = sqrt((xx-pos(lp,1)).^2 + (yy-pos(lp,2)).^2 + (zz-pos(lp,3)).^2);
    			rr(find(rr<radius)) = 10^9; % make it insignificant

			% Here's the change in frequency calculation:
			% cos(theta) at each location (theta is the angle
			% between position vector and Bo field):
    			costh = (zz-pos(lp,3))./rr;
			df = df + C*((radius./rr).^3).*(3*costh.^2 - 1);
		end
		% compute the mean of the signal equation for the voxel:
		sig(trial) = mean(exp(i*df(:)*2*pi*0.02)); 

	end
	msig(simn) = mean(abs(sig));
	sdsig(simn) = std(abs(sig));
	mph(simn) = mean(angle(sig));
	sdph(simn) = std(angle(sig));
	nb(simn) = nbub;
	br(simn) = radius;
	vf(simn) = vfrac;
	r2(simn) = -20/log(msig(simn));

	disp(sprintf('num = %d, rad = %f microns, vfrac = %g',nbub,radius,vfrac))
	disp(sprintf('Signal (TE = 20ms) mag = %f, phase = %f rad, T2p = %f',msig(simn),mph(simn),r2(simn)))

%end
