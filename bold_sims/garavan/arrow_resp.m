load times

% switch from scan units back to msec.
t_onset = onsets;
t_duration = lengths;

input=zeros(200 ,1);
for i=1:size(t_onset,1)
   input(floor(t_onset(i)):floor(t_onset(i) + t_duration(i)))=1;
end

delta=2500;
tau = 1250;
plot(input)
pause
gam=zeros(20000,1);
dummy=gam;
for t=1:size(gam,1)
   gam(t)=((t-delta)./tau).^2.*exp(-(t-delta)./tau).*((t-delta) > 0);
end
plot(gam);
pause

resp=conv(input,gam);
plot(resp)

pause

resp2=zeros(200,1);
for i=1:200
   resp2(i) = resp(i*2000);
end

plot(resp2);