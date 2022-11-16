u2 = [0 1 0]
b = [ 1 2 3]'
s = [10 20 30]



% this is the operation:


a = -j  * s  .* (u2  *  exp(-j*b*s))

% same thing without exp stuff
% a = u2 * s *  (b*s)

% we want the result to be:
-j*[10*exp(-j*2*10)    20*exp(-j*2*20)   30*exp(-j*2*30) ]


A = [ 1 2 ; 3 4]
B = [ 5 6 ; 7 8]
C = [9 10; 11 12]

A*(B.*C)

%%%

f = exp(-j*[1 2 3; 4 5 6; 7 8 9]) ;
b = [6 12 13]' * 1e-3
s = [4 5 6]
m = [7 8 9]'


f = exp(-j* rand(10,10));
b = rand(10,1)*1e-3;
s = sin(linspace(0,pi/5,10));
m = randn(10,1);


% forward model
y = (f.*exp(-j*b*s) )* m

% solving for b, by inverting things
bs = -angle ( (y*pinv(m)) ./ f)
 
bb = mean(bs,2)./s'

plot(b,bb,'o'); axis([-1 1 -1 1]) 
