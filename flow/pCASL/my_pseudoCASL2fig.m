%for fa=10:4:34
for fa=25    
npulses=1;

T1=1660;
T2=80;
%T1=1e6;
%T2=1e6;

dt=2e-3;  %dt in ms
endtime = 1.84;                 % ms
time = [0:dt:endtime-dt]';       % ms

gambar = 42570;               % gamma/2pi in kHz/T
gam=gambar*2*pi;

pw_rf1=.5;  %ms
mpgb1=hanning(round(pw_rf1/dt));

if(0)  % write the hanning pulse out to a binary file
b1pulse=mpgb1;
mymax=max(b1pulse);
b1pulse=floor((b1pulse./mymax.*32764)/2)*2;  %scale to full integer range
b1pulse(end)=b1pulse(end)+1;  %make the last number in the file odd
absname=('myhanning.bin');
absfid = fopen(absname ,'wb','b') ;  %must have the 'b' for big-endian byte ordering
fwrite(absfid,b1pulse,'short');  %must save as 16 bit integers to work with GE
fclose(absfid);
end


a_rf90=(fa/180*pi)/(trapz(abs(mpgb1))*dt*gam);
npts = length(time);          % number of time points for simulation
m = zeros([npts 3 ]);         % place holders
M = zeros(32,64);

%adavals = [0.02:0.02:.3];
%adavals=[0.12:.04:.28];
adavals=0.2;
%adavals=[.01:.01:.05];
%adavals=[.25:.05:.4];
%velocities = [6:2:60]; %cm/s
velocities=20;
for ada_index = 1:length(adavals)

clear beff gz 
beff = zeros([npts 3 ]);

ada=adavals(ada_index);  %fraction residual gradient moment

a_gz1 = 0.6/10000; % in T/cm  
a_gz2 = 6*0.6/10000; % in T/cm  (calculated to give zero moment)

a_gz2 = (1-ada)*a_gz2;

aramp=linspace(0,a_gz1,round(0.10/dt)).';
dramp=linspace(a_gz1,0,round(0.10/dt)).';
aramp2=linspace(0,-a_gz2,round(0.10/dt)).';
dramp2=linspace(-a_gz2,0,round(0.10/dt)).';

plateau1=a_gz1*ones(round((pw_rf1)/dt),1);
plateau2=a_gz1*ones(round((pw_rf1)/dt),1);
gap=zeros(round(0.02/dt),1);



gz =  [aramp; plateau1; dramp; aramp2; dramp2;  gap; aramp; plateau2; dramp; aramp2; dramp2; gap];
%gz =  [aramp; plateau1; dramp; gap; aramp; plateau2; dramp];
gz(end:length(time))=0;

zpos=10;

beff(round(0.1/dt):round(0.1/dt)+length(mpgb1)-1) = real(mpgb1)'.*a_rf90; 




for vel_index = 1:length(velocities)
clear Mz_ss* Mxy_ss* m m_ends
velocity = velocities(vel_index); %cm / s    


%control case
beff(round((0.1+pw_rf1+0.1+0.1+0.1+0.02+0.1)/dt):round((0.1+pw_rf1+0.1+0.1+0.1+0.02+0.1)/dt)+length(mpgb1)-1) = -real(mpgb1)'.*a_rf90; 
zstep = endtime/1000*velocity;
zposvals= -4:zstep:4;
m0 = [0 0 1];  % initial magnetization
m_ends(1,:)=m0;
for nn=1:length(zposvals)
zpos=zposvals(nn);   
beff(:,3) = gz.*zpos;  % + gx.*xpos + gy.*ypos;    
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

mfull=[];
%tag case
beff(round((0.1+pw_rf1+0.1+0.1+0.1+0.02+0.1)/dt):round((0.1+pw_rf1+0.1+0.1+0.1+0.02+0.1)/dt)+length(mpgb1)-1) = real(mpgb1)'.*a_rf90; 
m0 = [0 0 1];  % initial magnetization
m_ends(1,:)=m0;
for nn=1:length(zposvals)
zpos=zposvals(nn);   
beff(:,3) = gz.*zpos;  % + gx.*xpos + gy.*ypos;    
if(nn==1)
    m = blochsim(m0,beff,T1,T2,dt,npts);
else
    m = blochsim(m_ends(nn-1,:),beff,T1,T2,dt,npts);
end
m_ends(nn,:)=m(end,:);
mfull = [mfull; m];
Mz_ss_tag(nn)=m(end,3);
Mxy_ss_tag(nn)=m(end,1)+i*m(end,2);
end

%clear m m_ends
%figure,plot(zposvals*a_gz1*128.8e6,Mz_ss)
figure(1),
subplot(121),plot(zposvals,Mz_ss_control,zposvals,abs(Mxy_ss_control)),title('Control Case')
subplot(122),plot(zposvals,Mz_ss_tag,zposvals,abs(Mxy_ss_tag)),title('Tag Case')

% 
% for i=1:200; quiver3(0,0,0,m_ends(i,1),m_ends(i,2),m_ends(i,3)),axis([-1 1 -1 1 -1 1])
% ; hold off; drawnow; myMovie(i)=getframe; end


tag_efficiency = (1-Mz_ss_tag(end))/(2*exp(-(zposvals(end)/velocity)/(T1/1000)));  %correct for T1 recovery
control_efficiency = (1-Mz_ss_control(end))/(2*exp(-(zposvals(end)/velocity)/(T1/1000)));  %correct for T1 recovery
efficiency(ada_index,vel_index) = tag_efficiency-control_efficiency


end

end

%save pCASL1 adavals velocities efficiency
%

%eval(sprintf('save pCASLfa%d adavals velocities efficiency',fa));

if(0)
figure,
surf(velocities,adavals,efficiency),colormap jet
xlabel('Velocity (cm/s)')
ylabel('Fractional Residual Gradient Moment')
zlabel('Tagging Efficiency')
end


end
%figure,plot(adavals,effi



if(0)
    figure

    zvec1=linspace(zposvals(1),zposvals(end),size(mfull,1));
    for n=1:50:size(mfull,1)
        subplot(211)
        plot3([0 mfull(n,1)],[0 mfull(n,2)], [0 mfull(n,3)],'b','LineWidth',2),axis([-1 1 -1 1 -1 1])
        grid on
        view(90,0)
        hold on
        plot3([0 1],[0 0],[0 0],'k','LineWidth',1)
        plot3([0 0],[0 1],[0 0],'k','LineWidth',1)
        plot3([0 0],[0 0],[0 1],'k','LineWidth',1)
        hold off
        subplot(212)
        plot(zvec1,mfull(:,3),'b',zvec1(n),mfull(n,3),'rx','LineWidth',2)
        grid on
        drawnow
    end
end
RF_full=kron(ones(1,length(zposvals)).',beff);
RF_full=RF_full(:,1);


figure

zvec1=linspace(zposvals(1),zposvals(end),size(mfull,1));
for n=round(size(mfull,1)*7/16:10:size(mfull,1)*9/16)
    subplot(311)
    plot3([0 mfull(n,1)],[0 mfull(n,2)], [0 mfull(n,3)],'b','LineWidth',2),axis([-1 1 -1 1 -1 1])
    grid on
    view(90,0)
    hold on
    plot3([0 1],[0 0],[0 0],'k','LineWidth',1)
    plot3([0 0],[0 1],[0 0],'k','LineWidth',1)
    plot3([0 0],[0 0],[0 1],'k','LineWidth',1)
    hold off
    subplot(312)
    plot(zvec1,mfull(:,3),'b',zvec1(n),mfull(n,3),'rx','LineWidth',2)
    subplot(313)
    plot(zvec1(n-1000:n+1000),RF_full(n-1000:n+1000),'k',zvec1(n),RF_full(n),'rx','LineWidth',2)

    grid on
    drawnow
end





if(0)
M0=[0;0;1]
fa=25/180*pi;
R = [1 0 0; 
    0 cos(fa) sin(fa); 
    0 -sin(fa) cos(fa)]
ada=0.15;  %fractional moment
zloc=0.2; %cm
pa=ada*2*pi*4257*0.6*zloc*0.5e-3
Rz = [cos(pa) sin(pa) 0;
      -sin(pa) cos(pa) 0;
      0 0 1]

M1=R*M0;
M2=Rz*M1;
M3=R*M2;
M4=Rz*M3;
M5=R*M4;
M6=Rz*M5;



plot3([0 M2(1)],[0 M2(2)], [0 M2(3)],'g'),axis([-1 1 -1 1 -1 1])
hold on
plot3([0 M1(1)],[0 M1(2)], [0 M1(3)],'b'),axis([-1 1 -1 1 -1 1])
plot3([0 M3(1)],[0 M3(2)], [0 M3(3)],'r'),axis([-1 1 -1 1 -1 1])
plot3([0 M4(1)],[0 M4(2)], [0 M4(3)],'c'),axis([-1 1 -1 1 -1 1])
plot3([0 M5(1)],[0 M5(2)], [0 M5(3)],'m'),axis([-1 1 -1 1 -1 1])
plot3([0 M6(1)],[0 M6(2)], [0 M6(3)],'y'),axis([-1 1 -1 1 -1 1])
plot3([0 1],[0 0],[0 0],'k','LineWidth',2)
plot3([0 0],[0 1],[0 0],'k','LineWidth',2)
plot3([0 0],[0 0],[0 1],'k','LineWidth',2)
grid on
view(90,48)
end
