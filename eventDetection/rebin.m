function y = rebin(x, N)
% function y = rebin(x, N)
%
% x = input time series (binary)
% N = the number of neighbors on each side
% y = output time course
%---------
%   N = 2*N+1;
%   T = length(x);
% 	y = zeros(1,T/N)
% 	for c=1: N
% 		y = y + x(c:N:end);
% 	end
% 	y = y/N;
%
% return
%


	
N = 2*N+1;
T = length(x);
y = zeros(floor(T/N),1);
for c= 1:N
	y = y + x(c:N:N*length(y));
end
y = y/N;

return
	
