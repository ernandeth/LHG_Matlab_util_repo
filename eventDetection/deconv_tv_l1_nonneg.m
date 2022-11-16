function [xhat,A] = deconv_tv_l1_nonneg(y,H,h1,h2,x0)
%
% Usage: xhat = deconv_tv_l1(y,H,h1,h2)
%
% y is the data
% H is the HRF convolution matrix
% h1 is the L1 penalty weight
% h2 is the TV penalty weight
%
% estimates x from the cost function
% J(x) = ||Hx-y||^2 + h1*l1(x) + h2*TV(x) + alpha*N(x) 
%
% TV is the total variation penalty, l1 is the L1-norm penalty.
%
% last edits on 4/10/09 by Luis  (included the negativity penalty)
%
alpha = 100;
T=length(y);
rho=1e-8;
if(nargin==4)
    x0=y; % Initializer does not seem to matter
end
tol = 1e-5*T; %5

% save some time by doing these operations only once!
Hty = H'*y;
HtH = H'*H; 
D = diag(ones(T,1)) + diag(-ones(T-1,1),1); 
D=D(1:T-1,:);

for i=1:1000
	

	% total variation penalty: 
	delta_xt = diff(x0);
    w1 = h2/2 ./ sqrt(delta_xt.^2 + rho);
    W1 = diag(w1);

	% L1 penalty:
	w = h1/2./sqrt(x0.^2 + rho);
    
    W = diag(w);
    
    % negativity penalty:  
%     xtmp = x0;  xtmp(xtmp<0) = 0;
%     %w3 = h1/2 ./sqrt(xtmp.^2 + rho) ;
%     w3 = alpha ./sqrt(xtmp.^2 + rho) ;
    
    xtmp = x0;  xtmp(xtmp>0) = 0;
    w3 = alpha .* sqrt(xtmp.^2) ;
    
    W = diag(w3) + W;
    %%%
    
	% puts together the two penalties
	% speed this computation up: A = HtH + W + D'*W1*D;
    DtW1D = diag(w1(1:end)) + diag([0; w1(1:end-1)]) - diag(w1(1:end-1),1) - diag(w1(1:end-1),-1);
    DtW1D(T,T) = w1(end); 
	DtW1D(T-1,T) = -w1(end); 
	DtW1D(T,T-1) = -w1(end);

    % Least squares error :
	A = HtH + W + DtW1D;
	x1 = A \ (Hty);   % this the same as x1 = inv(A)*Hty , where inv(A) is computed by Gaussian elimination?
	
	% kludge: forcing solution to be positive
	% Bad idea - this didn't work:   x1 = abs(x1);  
	
	% res = y - H*x1;
    % if( abs(x1-x0)<1e-6 ) break; end
    if( norm(x1-x0) < tol ) break; end
    x0=x1;
end

fprintf('\r  Steps to convergence:  %d', i);
xhat=x1/max(x1);
xhat = xhat - mean(xhat);

return
