function G =dictionary_correlation_150507 (x)
% function G = dictionary_correlation_150507(x)
%
% this function  takes in a vector of timing parameters and uses it to
% create a dictionary of ASL signals.
% it returns the correlation among the dictionary entries
%
%
% it returns the norm of the correlation matrix of the dictionary
%

% duplicate the timing vector so that we can do teh control and the tag the
% same. 
%tvec = [x x];
%tvec(:) =0;
%tvec(1:2:end) = x;
%tvec(2:2:end) = x;

% no longer doing things in pairs.
tvec = x;
%maxTR = 0.7;

N = length(tvec)/2;

% %
% timing_parms.label_dur =  tvec(1 : N);
% timing_parms.PID = tvec(N+1 : 2*N) ;
% %timing_parms.t_adjust =   tvec(2*N+1 : 3*N );
% timing_parms.t_adjust =  maxTR-  tvec(1 : N) - tvec(N+1 : 2*N) ;
% timing_parms.order = 1;

% timing_parms.PID = tvec(N+1 : 2*N) ;
% timing_parms.t_adjust =   tvec(2*N+1 : 3*N );

timing_parms.label_dur =  tvec(1 : N) ;
timing_parms.PID = 0.05 * ones(1,N) ;
timing_parms.t_adjust =0.05 * ones(1,N) ;
timing_parms.isLabel = tvec(N+1 : 2*N)  ;
timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
 
%
%timing_parms
[dict, parms] = gen_dictionary_150507(timing_parms);

%plot(dict')
%drawnow

xc = corrcoef(dict');
% I only care about half of the matrix
for n=1:size(xc,1)
    xc(n, 1:n) = 0;
end

% penalize entries greater than 0.9 correlation
xc(xc>0.9) = 1.5*xc(xc>0.9);

xc = (xc(:));
G = mean(xc);
%G = norm(xc);
return
    
    
