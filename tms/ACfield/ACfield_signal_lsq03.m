function result = ACfield_signal_lsq03(Bn, parms, yn)
% function ACfield_signal_lsq03(Bi, parms, yn)
%
% This function is used to solve for B univariate equation.
% It computes the difference between the RHS and LHS of an equation
% the parent function tries to minimize the difference.
%
% LHS = sum(yn * Fn .* S .* sin(Bn*S) .* M0);
% 
% RHS = 0;
% for kp=1:length(M0)
%     
%     RHS = RHS + sum(...
%         (Fn * Fn(kp)) .* ...
%         sin(Bn.*(S - S(kp) )) .* ...
%         (S-S(kp)) .* M0 * M0(kp) );
% end
% 
%  if called with two parameters only, this function generates an MR sample
%
% if a third input parameter is added (yn), we assume that it's been
% called by a fitting routing, like lsqnonlin.  In that case, the function returns the RSS
% difference between those yn and a signal generated from input parameters.
%
% all the input parameters are stuck into the structure 'parms'
%


%global parms
Fn = parms.Fn ;    % nth row of the FFT matrix
M0 = parms.M0;    % the object
S =  parms.ACfun;  % integral of B field waveform (in time) : phase contribution in rads.

Fn = Fn(:);
M0 = M0(:);
S = S(:);

% old derivation was wrong: yn = (Fn .* exp(-j*Bn*S) ) * M0;


LHS = sum(yn * Fn .* S .* sin(Bn*S) .* M0);

RHS = 0;

for k=1:length(M0)
        RHS = RHS + sum(...
            (Fn * Fn(k)) .* ...
            sin(Bn*(S(k) - S )) .* ... 
            (S-S(k)) .* (M0 * M0(k)) );
end
RHS = 0.5*RHS;

%if nargin>2
result = norm(LHS -  RHS) ;
%else
%    result = (Fn.*exp(-j*Bn*S))*M0 ;
%end


return

