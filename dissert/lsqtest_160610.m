function lsqtest

p.x1 = ones(100,1);
p.x2 = sin(linspace(0,5*pi,100))';


data = stuff([1 2], p);

data = data + 0.5*randn(size(data));
hold off, plot(data , 'o')

myfun = @(guess)stuff(guess,p, data);

opts = optimset('lsqnonlin');

% NOTE that the data does NOT get entered here
estimates = lsqnonlin(myfun , [0 0], [-10 -10], [10 10], opts)

thefit = stuff(estimates, p);

hold on, plot(thefit)
return




function result = stuff(c, p, data)


x1 = p.x1;
x2 = p.x2;

y = c(1)*x1 + c(2)*x2;

if nargin==3
    result = norm(data - y);
else
    result = y;
end

return
