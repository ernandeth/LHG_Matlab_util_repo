function result = ps_lsq(guess, TE, data)
% result = ps_lsq(guess,  TE, data)
%
% defines the relaxation in an Progressive Saturation experiment.
% Note: uses the magnitude images (no phase)
% Mo = guess(1);
% T2 = guess(2);
%   result = Mo *exp(- TE / T2) ;

    Mo = guess(1);
    T2 = guess(2);
    
    t = TE;
    %t = data;
    
   M_guess = Mo *exp(-TE / T2) ;
    
    if nargin==2
        result = M_guess;
    else
        M = data;
        result = (M - M_guess);
    end
               
return
