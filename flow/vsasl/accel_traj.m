seg1 = 0.1*sin(linspace(0,2*pi, 50));
seg1 = ones(1,50);
seg1(end/2+1:end) = -1;
seg2 = zeros(size(seg1));
module = [seg1 -seg1 seg2 seg2];

rf1 = sinc(linspace(-2*pi, 2*pi, 100));
rf = [seg2 seg2 rf1 ];

g = repmat(module, 1,5);
rf = repmat(rf, 1,5);

t = linspace(0,10,length(g));
kr = cumsum(g);
kv = cumsum(t.*g);
ka = cumsum(t.*t.*g);

zmask = ones(size(g));
zmask(abs(g) > eps) = 0;

% k trajectorie
subplot(421)
plot(t,rf,'r'), grid on
title('RF pulses')
hold on
area(t,g)
title('gradient')
hold off

subplot(423)
plot(t,kr,'k'), grid on
title('k_r')
hold on
plot(t,rf,'r')
hold off

subplot(425)
plot(t,kv,'k'), grid on
title('k_v')
hold on
plot(t,rf,'r')
hold off

subplot(427)
plot(t,ka,'k'), grid on
title('k_a')
hold on
plot(t,rf,'r')
hold off

% Zoomed in version:  clip out the trajectory during gradient readout
subplot(422)
area(t,g), grid on
title('gradient')
hold on
plot(t,rf,'r')
hold off

subplot(424)
plot(t, kr .* zmask ,'k'), grid on
title('k_r')
hold on
plot(t,rf,'r')
hold off

subplot(426)
plot(t, kv.*zmask,'k'), grid on
title('k_v')
hold on
plot(t,rf,'r')
hold off

subplot(428)
plot(t, ka.*zmask,'k'), grid on
title('k_a')
hold on
plot(t,rf,'r')
hold off




accel = (50-4)/1.5
