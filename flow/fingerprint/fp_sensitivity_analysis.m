% fp_sensitivity script
names = { 'mtis0',     'f',  'cbva' ,  'bat',  'bat2', 'kfor',  'r1tis',  'flip',  'Disp'};



sig1 = gen_signals_180320(parms,  timing_parms, 0,0) ;
% first look at the sensitivity of the signal
% norm of each partial - df/f0
delts = calc_partials(timing_parms, parms, 1);

% we don't care about M0 and bat2
delts = delts(2:end, :);
delts(4, :) = [];
sensitivity = sum(delts.^2, 2) / size(delts,2);


% Now look at the CRLB
% calculates the actual partial derivatives:  df/dparm
partials = calc_partials(timing_parms, parms, 0);
% we don't care about M0 and bat2 and Disp
names(1) = [];
names(4) = [];
names(end) = [];

partials(1, :) = [];
partials(4, :) = [];
partials(end, :) = [];

p=struct2arr(parms);
p(1) = [];
p(4) = [];
p(end) = [];

% estimate of the noise
sigma2 = 0.01;

% Fisher information matrix
FImat = (1/sigma2) * partials * partials' ;

% CLRB is the inverse of the FI matrix
CRLB = inv(FImat);


%names

 fprintf('\nCRLB: \n');
 fprintf('%3g\t', 100*sqrt(diag(CRLB))' ./ p' );
%  fprintf('\nsensitivity: \n');
%  fprintf('%.3g\t', sensitivity');
 


 figure
 subplot(311)
 plot(timing_parms.t_tag)
 title('Tagging Duration schedule')
 subplot(312)
 plot(sig1)
 title('Signal')
 subplot(313)
 plot(delts(1:4,:)'); legend(names(1:4))
 axis([0 size(delts,2) -5e-3 5e-3])
 title('Signal Change WRT 10% change in parameter...')
 
