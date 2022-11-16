function gout = Ggen(shape,G,r,f,gap,dt);

% G in G/cm
% r raise time, in ms
% flat time, in ms
% gap before, after and between gradients, in ms
% dt in us
% Jia Guo, 2011@ucsd

dt = dt*1e-6; % convert temporal resolution from us to s
% convert ms to s
r = r/1000;
f = f/1000;
gap = gap/1000; % pad with 0 at the end
npulse =  numel(shape);
nr = round(r/dt);
nf = round(f/dt);
ngap = round(gap/dt);
gout = zeros(1,ngap(1));
for i=1:npulse
    tmp = [];
    switch shape(i)
        case 't' %trapezoid
            tmp = [G(i)*[1:nr(i)]/nr(i), G(i)*ones(1,nf(i)), G(i)*[nr(i):-1:1]/nr(i)];
        case 's' % sin
            tmp = [G(i)*sin([1:ceil(nf(i)/2)]/ceil(nf(i)/2)*pi/2), G(i)*sin([ceil(nf(i)/2):-1:1]/ceil(nf(i)/2)*pi/2)];
        case 'r' % rect
            tmp = [G(i)*ones(1,nf(i))];
        case 'c'
            tmp = G(i)*(1-cos([1:ceil(nf(i))]/ceil(nf(i))*pi/2));
            tmp = [tmp, tmp(end:-1:1)];
    end
    gout = [gout tmp zeros(1,ngap(i+1))];
end