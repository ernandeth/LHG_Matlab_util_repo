% comapre HRF parameters

tau1=17;
tau2=42;
T = 20;
p = T-10;
TP = 20;


h = inline(' (t).^(tau1).*exp(-t)/gamma(tau1) - 0.1*(t).^(tau2).*exp(-t)/gamma(tau2)',...
    'tau1','tau2','t');
t = [0:2:2*T-2]';


for h_err=linspace(-2,2,3)
hh = h(tau1+h_err,tau2, t+7);
hh = hh - mean(hh);
hh = hh/norm(hh);
find(hh==max(hh))
plot(hh, 'k'); hold on

end

title('Variabilty of HRF')
xlabel('Time (s)')
ylabel('Intensity (a.u.)')
dofontsize(16)
fatlines