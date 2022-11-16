function result = ACfield_signal_lsq(ACfield, parms, data)
% function ACfield_signal_lsq(ACfield, parms, data)
%
% if called with two parameters only, this function generates an MR signal
% based on the signal equation with an extra term that accounts for the
% extra phase induced by an external oscillating magnetic field
%
% if a third input parameter is added (data), we assume that it's been
% called by a fitting routing, like lsqnonlin.  In that case, the function returns the RSS
% difference between those data and a signal generated from input parameters.
%
% all the input parameters are stuck into the structure 'parms'
%
% M0_kphase = parms.M0kphase;  % this should be an r x t matrix
% ACphase = parms.ACphase;     % this should be an r x t matrix
%
% This is how those were computed:
%
% M0_phase(r,:) = M0(r) * exp(i * kxy * r) * exp (i * gammabar * Boff(r));
% ACphase(:,t) = gammabar * cumsum( sin(2*pi*freq*t + initial_phase*pi/180) );
%

%global parms

M0_kphase = parms.M0phase;  % this should be an r x t matrix
ACphase = parms.ACphase;     % this should be an r x t matrix
inds = parms.inds;

DEBUG=0;
if DEBUG
    xdim = sqrt(size(M0_kphase,1));
    imagesc(reshape(ACfield,xdim, xdim)); drawnow
end

% original code from Vivek:
% signal...
% 	= sum( (parms.m(ind) * ones(1,length(parms.phase_vec)) ...
% 	.* exp(-1i .* 2 .* pi .* (parms.yy_mat(ind) * parms.kro + parms.xx_mat(ind) * parms.kpe))...
% 	.* exp(1i .* parms.gambar .* Bfield_est * parms.phase_vec')), 1 );

% New version:  a bunch of stuff is pre-computed
% calculate the phase that each pixel gets because of the magnitude of the AC field
% at that pixel:
%
% ACphase_r = kron(ACfield(:), ACphase);
%
% now include that AC field phase into the contribution of that pixel to the signal
% and compute the integral
%
% signal = double(sum(M0_kphase .* exp( i * ACphase_r ) ,1 ));
%

% To save time, we only compute the pixels there is enough spin density
signal = 0;
for r=inds'
    signal = ...
        signal +  ...
        M0_kphase(r,:).* exp(i* ACfield(r)*ACphase) ;
end
signal=double(signal);


if nargin>2
    result = sum(abs(data -  signal)) ; % + (parms.beta * norm(gx(:)+gy(:)));
else
    result = signal;
end



return

