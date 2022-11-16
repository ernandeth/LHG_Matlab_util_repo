fa=20;
npulses=1;

%T1=1660;
%T2=140;
T1=1e6;
T2=1e6;

dt=2e-3;  %dt in ms
%endtime = 1.84-2*.45;                 % ms
endtime = 1.84;
time = [0:dt:endtime-dt]';       % ms

gambar = 42570;               % gamma/2pi in kHz/T
gam=gambar*2*pi;

%pw_rf1=0.05;  %ms
pw_rf1 =0.5; %ms
mpgb1 =hanning(round(pw_rf1/dt));
a_rf90 =(fa/180*pi)/(trapz(abs(mpgb1))*dt*gam);
npts = length(time);          % number of time points for simulation
m = zeros([npts 3 ]);         % place holders
M = zeros(32,64);
beff = zeros([npts 3 ]);

a_gz1 = 0.6/10000; % in T/cm  
a_gz2c = 3.6/10000; % in T/cm  
%a_gz2c = 0.9/10000; % in T/cm  
ada=0.11;  %fraction residual gradient moment
a_gz2 = (1-ada)*a_gz2c;
a_gz2c = a_gz2;

aramp = linspace(0,a_gz1,round(0.10/dt)).';
dramp = linspace(a_gz1,0,round(0.10/dt)).';
aramp2 = linspace(0,-a_gz2,round(0.10/dt)).';
dramp2 = linspace(-a_gz2,0,round(0.10/dt)).';
aramp2c = linspace(0,-a_gz2c,round(0.10/dt)).';
dramp2c = linspace(-a_gz2c,0,round(0.10/dt)).';

plateau1 = a_gz1*ones(round((pw_rf1)/dt),1);
plateau2 = a_gz1*ones(round((pw_rf1)/dt),1);
gap = zeros(round(0.02/dt),1);

gz =  [aramp; plateau1; dramp; aramp2; dramp2;  gap; aramp; plateau2; dramp; aramp2; dramp2; gap];
gzc =  [aramp; plateau1; dramp; aramp2c; dramp2c;  gap; aramp; plateau2; dramp; aramp2c; dramp2c; gap];
%gz =  [aramp; plateau1; dramp; gap; aramp; plateau2; dramp];
gz(end:length(time))=0;

zpos=10;

beff(round(0.1/dt):round(0.1/dt)+length(mpgb1)-1) = real(mpgb1)'.*a_rf90; 

ang_tmp=trapz([aramp; plateau1; dramp])*dt*gam*180/pi*ada


RFcenter=1; %cm
B1offset=RFcenter*2.554/gambar; %Tesla  (assumes 0.6G/cm gradient, thus 2.554 kHz / cm)



velocity = 20; %cm / s

zstart=-15;
zend=15;


%control case
beff(round((0.1+pw_rf1+0.1+0.1+0.1+0.02+0.1)/dt):round((0.1+pw_rf1+0.1+0.1+0.1+0.02+0.1)/dt)+length(mpgb1)-1) = -real(mpgb1)'.*a_rf90; 
zstep = endtime/1000*velocity;
zposvals= zstart:zstep:zend;
m0 = [0 0 1];  % initial magnetization
m_ends(1,:)=m0;
for nn=1:length(zposvals)
zpos=zposvals(nn);   
beff(:,3) = gzc.*zpos+(beff(:,1)>0)*B1offset;  % + gx.*xpos + gy.*ypos;    
if(nn==1)
    m = blochsim(m0,beff,T1,T2,dt,npts);
else
    m = blochsim(m_ends(nn-1,:),beff,T1,T2,dt,npts);
end
m_ends(nn,:)=m(end,:);
Mz_ss_control(nn)=m(end,3);
Mxy_ss_control(nn)=m(end,1)+i*m(end,2);
end

clear m m_ends

%tag case
beff(round((0.1+pw_rf1+0.1+0.1+0.1+0.02+0.1)/dt):round((0.1+pw_rf1+0.1+0.1+0.1+0.02+0.1)/dt)+length(mpgb1)-1) = real(mpgb1)'.*a_rf90; 
zstep = endtime/1000*velocity;
zposvals= zstart:zstep:zend;
m0 = [0 0 1];  % initial magnetization
m_ends(1,:)=m0;
for nn=1:length(zposvals)
zpos=zposvals(nn);   
beff(:,3) = gz.*zpos+(beff(:,1)>0)*B1offset;  % + gx.*xpos + gy.*ypos;    
if(nn==1)
    m = blochsim(m0,beff,T1,T2,dt,npts);
else
    m = blochsim(m_ends(nn-1,:),beff,T1,T2,dt,npts);
end
m_ends(nn,:)=m(end,:);
Mz_ss_tag(nn)=m(end,3);
Mxy_ss_tag(nn)=m(end,1)+i*m(end,2);
end


%figure,plot(zposvals*a_gz1*128.8e6,Mz_ss)
figure,
subplot(121),plot(zposvals,Mz_ss_control,zposvals,abs(Mxy_ss_control)),title('Control Case')
subplot(122),plot(zposvals,Mz_ss_tag,zposvals,abs(Mxy_ss_tag)),title('Tag Case')

% 
% for i=1:200; quiver3(0,0,0,m_ends(i,1),m_ends(i,2),m_ends(i,3)),axis([-1 1 -1 1 -1 1])
% ; hold off; drawnow; myMovie(i)=getframe; end


tag_efficiency = (1-Mz_ss_tag(end))/(2*exp(-(zposvals(end)/velocity)/(T1/1000)));  %correct for T1 recovery
control_efficiency = (1-Mz_ss_control(end))/(2*exp(-(zposvals(end)/velocity)/(T1/1000)));  %correct for T1 recovery
efficiency = tag_efficiency-control_efficiency
