function df = cdiff(x)
% function df = cdiff(x)
% central difference 
% output is same size as input.
% uses zero padding for this.
x = x(:);
x = [0;x;0];
df = (x(2:end-1)-x(1:end-2) + x(3:end) -x(2:end-1))/2;

return