function output=sinc(in)
% does sinc(x) function sin( pi .* x) ./ (pi .* x)
output = sin(pi.*(in+eps)) ./ (pi.*(in+eps));
