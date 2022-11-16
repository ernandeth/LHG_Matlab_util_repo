



function junk

x=3;
A=4;

data = exp(-i*A*x);

mypol(3, A, data)

plot(exp(-i*A*[-3:0.1:3]))

fzero(@(xhat) mypol(xhat, A, data), 1)

return

function y = mypol(x, A, data)

model = exp(-i*A*x);
y = real(model) - real(data) - imag(model) + imag(data);

return
