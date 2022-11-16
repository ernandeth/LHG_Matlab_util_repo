% Generic min-max temporal interpolators

t = [0:4000].*5e-6;
tmax = max(t);
L = 10;
tau = tmax/L;
t2 = [0:L]*tau;
nn = 0*(randn(size(t2)) + i*randn(size(t2)));
df = 100;
[ll llt] = meshgrid([0:L]);
ggt2 = sinc((ll-llt)*tau*df).^2;
gb2 = sinc((([0:L]'*tau*ones(size(t)) - ones([L+1 1])*t)*df)).^2;
beta = 0; %.00001;
a2 = inv(ggt2 + beta*eye(L+1))*gb2;
figure(1)
subplot(411)
plot(real(a2'))
subplot(412)
plot(abs((real((ones(size(t2))+nn)*a2))-1))
subplot(413)
f1 = 50;
plot(abs((exp(i*2*pi*f1*t2)+nn)*a2 - exp(i*2*pi*f1*t)));
subplot(414)
f1 = 100;
plot(abs((exp(i*2*pi*f1*t2)+nn)*a2 - exp(i*2*pi*f1*t)));

if 0
figure(2) 
beta = 0.00001;
a2 = inv(ggt2 + beta*eye(L+1))*gb2;
subplot(411)
plot(real(a2'))
subplot(412)
plot(abs((real((ones(size(t2))+nn)*a2))-1))
subplot(413)
f1 = 50;
plot(abs((exp(i*2*pi*f1*t2)+nn)*a2 - exp(i*2*pi*f1*t)));
subplot(414)
f1 = 100;
plot(abs((exp(i*2*pi*f1*t2)+nn)*a2 - exp(i*2*pi*f1*t)));
return
end
% add the ^2 for triangle histogram
% it turns out the rect(200) has nearly the same interpolators as tri(100)
df = 200;

nb = 41;
f = linspace(-df/2,df/2,nb);
hist = ones(size(f));
hist = df/2-abs(f);

H = diag(hist)./sum(hist);
P = exp(i*2*pi*t2'*f).';
E = exp(i*2*pi*t'*f).';
figure(2)
a3 = inv(P'*(H)*P + beta*eye(L+1))*P'*H*E;
subplot(411)
plot(real(a3'))
subplot(412)
plot(abs(sum(real(a3))-1))
subplot(413)
f1 = 50;
plot(abs(exp(i*2*pi*f1*t2)*a3 - exp(i*2*pi*f1*t)));
subplot(414)
f1 = 100;
plot(abs(exp(i*2*pi*f1*t2)*a3 - exp(i*2*pi*f1*t)));

return

% Brad's MIN-MAX
for lp1 = 0:L
    for lp2 = 0:L
        ggt(lp1+1,lp2+1) = sum(hist.*exp(-i*2*pi*f*tau*(lp1-lp2)));
    end
end
ggt = ggt./sum(hist);

gb = zeros([L+1 length(t)]);
for lp1 = 0:L
    gb(lp1+1,:) = sum((hist'*ones(size(t))).*exp(-i*2*pi*f'*(t-tau*lp1)),1);
end
gb = gb./sum(hist);

a = inv(ggt)*gb;

