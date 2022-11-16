function Gmap = makeGmap(S, psi)
% function Gmap = makeGmap(S, psi)
%
% this computes the G factor map from sensitivitiy maps and noise
% correlation matrix from Pruessmann:  Magnetic Resonance in Medicine 42:952?962 (1999)
% 
%      Gmap(:) = diag(sqrt( inv(SH * psi_inv * S) .* (SH * psi_inv * S)));

Ncoils = size(S,1);
Npix = size(S,2);

psi_inv = pinv(psi);
%psi_inv=eye(8);


SH = S';
%SH = transpose(conj(S));

Gmap = sqrt(diag( pinv(SH * psi_inv * S)) .* diag(SH * psi_inv * S)) ;
Gmap = real(Gmap);
%Gmap = reshape(Gmap, sqrt(Npix), sqrt(Npix));
%figure(3); imagesc(abs(Gmap));


return