function result = BOLD_lsq(guess, t, B )
% function result = BOLD_lsq(guess, t [,B])
%
% defines the gamma variate BOLD response function for 
% fitting purposes.
%

   

    delta = guess(1);
    Tau= 	guess(2);
	Amp =	guess(3);
     
   
   B_guess=zeros(size(t));

   B_guess= Amp* make_hrf(delta, Tau, max(size(t)) );
    
    if nargin==3
        result = (B - B_guess) ;
    else
        result=B_guess;        
    end
    

return
