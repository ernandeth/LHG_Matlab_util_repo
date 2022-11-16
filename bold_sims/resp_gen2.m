function resp2 = resp_gen2(TR,DURATION, ISI)
% function resp2 = resp_gen2(TR, DURATION, [,ISI})
%
% Luis Hernandez
% University of Michigan
%
% This program generates a simulated BOLD response
% to a set of EVENTS that occur every 50 secs, 
% or at a given ISI, if specified.  starting at
% a particular delay (hence the input variable 'delay')
%
% The total imaging time is DURATION sec.
%
% The number of points in one second is specified by the 
% variable 'points'.( not related to TR )
%

points = 10; % number of points in one second
delay = 0;

if nargin == 1
    ISI=50;
end

stim=zeros(DURATION*points,1);
stim(1+ delay*points :  ISI*points :  DURATION*points) =1;


%hrf = make_hrf(2.5*points,1.25*points, 20*points);
hrf = spm_hrf(1/points );

resp=conv(stim,hrf);
resp = resp(1:size(stim));

subplot 311, plot(hrf)
subplot 312, plot(stim)
subplot 313, plot(resp)


% resample the data to the acquisition rate 
TR = TR*points;

stim2 = stim(1:TR:size(stim));
hrf2 = hrf(1:TR:size(hrf));
resp2 = resp(1:TR:size(resp));

%disp('press return to continue')
%pause

subplot 311, plot(hrf2)
subplot 312, plot(stim2)
subplot 313, plot(resp2)

return

%%%%%
function hrf = make_hrf(delta,tau, npoints)
% hrf = make_hrf(delta,tau, npoints)
%
% Creates a BOLD - Hemodynamic Response Function simulated by a 
% gamma-variate function
%
% 'npoints' is the number of points in the data array
% (very important that everything is in the same units)
% 
hrf=zeros(npoints,1);

	for t=1:size(hrf,1)
   		hrf(t)=((t-delta)./tau).^2.*exp(-(t-delta)./tau).*((t-delta) > 0);
   end  
      
return

%%%%
function result = make_squares(onsets, lengths)
% function result = make_squares(onsets, lengths)
%
% creates a square wave function 
% 'onsets' is a vector of onset times
% 'lengths' is a vector of the durations of the squares
%

	t_onset = onsets;
	t_duration = lengths;

	result=zeros(200 ,1);
	for i=1:size(t_onset,1)
   		result(floor(t_onset(i)):floor(t_onset(i) + t_duration(i)))=1;
	end
   
return
