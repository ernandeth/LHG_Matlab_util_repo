
fixtime= [7 9 9 9 9 9 10 10 10 10 10 10 5 5 5 5 5 5 7 7 7 7 7 7 6 6 6 6 6 6 8 8 8 8 8 8 ]+1;
onsets = cumsum(fixtime);

y=load('../events_sphere_tdata.dat');


[hrf v] = hrf_deconv(y,onsets, 30);

hrf = hrf-mean(hrf);
hrf = hrf/max(hrf);

plot(hrf);

tau1=18; tau2=30;
H = HRF_mat(tau1, tau2, 30);

hold on
plot(H(:,1),'g');



%%\
tc = load('../BlockEvent01_sphere_tdata.dat');
tc = mydetrend(tc);
TR=1
duration = 255 * TR;
fixtime= [9 9 9 10 10 10 10 10 10 7 7 7 6 6 6 6 6 6 8 8 8 0];
acttime = [0 10 10 10 1 1 1 1 1 1 5 5 5 1 1 1 1 1 1 10 10 10];

onsets = fixtime + acttime;
onsets = cumsum(onsets);

z=zeros(size(tc));
z(onsets)=10;

onsets = onsets([4:9 15:18]);
z2=zeros(size(tc));
z2(onsets)=10;
hrf=0;
for c=1:length(onsets)
    tmp = tc(onsets(c):onsets(c)+15);
    %tmp = tmp - mean(tmp);
    hrf = hrf + tmp;
end
hrf = hrf - mean(hrf);
hrf = hrf/max(hrf);

plot(hrf);

tau1=12; tau2=30;
H = HRF_mat(tau1, tau2, 15);

hold on
plot(H(:,1),'g');
