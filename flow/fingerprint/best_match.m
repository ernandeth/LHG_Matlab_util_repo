function [best r]= best_match(dict, y)
% function [best, r] = best_match(dict, y)
%
% search dictionary for best match
% the match is determined by inner product
% ---OR ----
% by the correlation coefficient
% (check the code)
%

%innerp = zeros(size(dict,2),1);
% for d=1:size(dict,2)
%     %innerp(d) = corr(y, dict(:,d)); 
%    
%     % 2/20/15
%    % restore this:
%     innerp(d) = y' * dict(:,d); 
% end

% Get rid of the for loop !!!!
innerp = y' * (dict);


%{
% mean centering:
y = y.' - mean(y);

dict_temp = dict  - repmat(mean(dict,1) , [length(y) 1] );

% normalizing 
dict_temp_norm = sqrt(sum(dict_temp .* conj(dict_temp)));

% computing inner products of the normalized signals.
innerp = conj(y) * dict_temp ./ dict_temp_norm ./ norm(y);
%}


% find the best match in the reduced dictioanry (highest inner product)
%best = find((innerp) == r);


% 2/20/15
[r best] = max(abs(innerp));

return

%% ----------- testing code -----
phys_parms.f =         0.008;
phys_parms.mtis0 =     1 ;
phys_parms.cbva =      0.025 ;
phys_parms.transit=    0.45 ;
phys_parms.kfor =      0.3; % 1e-2 ;
phys_parms.r1tis =     1/1.7  ;
phys_parms.beta =      75*pi/180 ; % flip angle in radians
phys_parms.L = 1;
phys_parms.Disp =      13;

for Disp = [1:5:30]
    phys_parms.Disp = Disp
% Just a quick test to make sure everything works as expected
obs = gen_signals_140328(phys_parms, timing_parms, 0, 1);
obs = obs / norm(obs);
plot(obs,'k')
drawnow
end
axis([1 25 -0.07 0.1])

