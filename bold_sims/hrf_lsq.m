function result = hrf_lsq(guess, t, data)
% function result = hrf_lsq( guess_params, t [, data])
%
% Generates a gamma variate function for the hemodynamic response function.
% points is the number of points to be included in the function.
% (ie - the length of the response to be calculated)
% Normalized so that it peaks at 1
%
%
% delta 	- delay in the function
% tau		- decay/uptake constant
% amp 	- amplitude  function
% 
% 
% delta = guess(1);
% tau = guess(2);
% 
% delta2 = guess(4);
% tau2 = guess(5);
% amp2 = guess(6);
% 
% % the amplitude is assumed to be 1
% t = a time vector corresponding to the data samples)
% data 		- the data to be fit
%
% NOTE:  must be consistent in units!
	   
delta = guess(1);
tau = guess(2);

delta2 = guess(3);
tau2 = guess(4);
amp2 = guess(5);

hrf_guess =  ((t-delta)./tau).^2  .*  exp(-(t-delta)./tau)  .* ((t-delta) > 0) - ...
	 amp2*((t-delta2)./tau2).^2  .*  exp(-(t-delta2)./tau2)  .* ((t-delta2) > 0) ;
				
hrf_guess =  hrf_guess./(max(abs(hrf_guess)));


if nargin==3
	weights = ones(size(data)); % remove ?
	weights(round(end*0.5):end) = 0.5; % remove ?
	result = data - hrf_guess;
	result = result .* weights ;  % remove ?
else
	result = hrf_guess;
end

return
