
% script to generate VSI pulses
% Jia Guo@UCR, 2019

%% generate RF, G pulses
dt = 1;
gamma = 4257; %Hz/Gauss

% RF
rf = [];
pw180 = 0.5; % ms, 180 FA
pw90 = 0.248; % 90 FA
pwSFA = 0.4; %ms, small flip angle RF
pwcomp = pw180 + pw90*2;
n180 = round(pw180*1000/dt);
n90 = round(pw90*1000/dt);
nSFA = round(pwSFA*1000/dt);
pw180 = 0.5; % ms, 180 FA
        
% sinc modulation
nmod = 8;
step = pi/(nmod/2+1);
SFA = sinc(-pi+step:step:pi-step);      

% fix b1
b1 = pi/(2*pi*gamma)/(pw180*1e-3); % composite 180
b1sc = 0.5/b1/gamma/(pwSFA/1000)/sum(SFA); % for small flip angle pulses

% Gradient
rtime = 0.3; % ms, rise time in gradient
ftime = 0.3; % ms, flat time in gradient
vspad1 = 0.08; % ms, gap before gradient
vspad2 = 0.42; % ms, gap after gradient
gtime = 2*rtime+ftime; %ms
ng = round(1000*gtime/dt);

% composite pulse
P90a = []; % 90 FA, 0 phase
P90a(1,:) = ones(1,n90);
P90a(2,:) = zeros(1,n90);

P90b = []; % 90 FA, pi phase
P90b(1,:) = ones(1,n90);
P90b(2,:) = pi*ones(1,n90);

P180a = []; % 180 FA, pi/2 phase
P180a(1,:) = ones(1,n180);
P180a(2,:) = pi/2*ones(1,n180);

P180b = []; % 180 FA, -pi/2 phase
P180b(1,:) = ones(1,n180);
P180b(2,:) = -pi/2*ones(1,n180);

comp1 = [P90a,P180a,P90a];
comp2 = [P90b,P180b,P90b];

PSFA0 = [];
PSFA0(1,:) = b1sc*ones(1,nSFA);
PSFA0(2,:) = zeros(1,nSFA);
            
P90a = []; % 90 FA, 0 phase
P90a(1,:) = ones(1,n90);
P90a(2,:) = zeros(1,n90);

P90b = []; % 90 FA, pi phase
P90b(1,:) = ones(1,n90);
P90b(2,:) = pi*ones(1,n90);

P180a = []; % 180 FA, pi/2 phase
P180a(1,:) = ones(1,n180);
P180a(2,:) = pi/2*ones(1,n180);

P180b = []; % 180 FA, -pi/2 phase
P180b(1,:) = ones(1,n180);
P180b(2,:) = -pi/2*ones(1,n180);

comp1 = [P90a,P180a,P90a];
comp2 = [P90b,P180b,P90b];

PSFA0 = [];
PSFA0(1,:) = ones(1,nSFA);
PSFA0(2,:) = zeros(1,nSFA);

g0 = Ggen('tttt',[1 -1 1 -1],rtime*[1 1 1 1],ftime*[1 1 1 1],[pwSFA+vspad1, pwcomp+vspad1+vspad2, pwSFA+vspad1+vspad2, pwcomp+vspad1+vspad2, vspad2],dt);
g = [repmat(g0,[1,8]),zeros(1,nSFA)];

gtime_gap = 2*rtime+ftime + vspad1 + vspad2; %ms
ngtime_gap = round(1000*gtime_gap/dt);

rf = [  SFA(1)*PSFA0, zeros(2,ngtime_gap), comp1, zeros(2,2*ngtime_gap+nSFA), ...
    comp1, zeros(2,ngtime_gap), SFA(2)*PSFA0,... 
    zeros(2,ngtime_gap),... 
    comp2, zeros(2,2*ngtime_gap+nSFA),... 
    comp2, zeros(2,ngtime_gap), ...
    SFA(3)*PSFA0, ...
    zeros(2,ngtime_gap), ...
    comp2, ...
    zeros(2,2*ngtime_gap+nSFA),... 
    comp1, ...
    zeros(2,ngtime_gap), ...
    SFA(4)*PSFA0, ...
    zeros(2,ngtime_gap), ...
    comp1, ...
    zeros(2,2*ngtime_gap+nSFA), comp2, zeros(2,ngtime_gap), ...
    SFA(5)*PSFA0, zeros(2,ngtime_gap), comp2, zeros(2,2*ngtime_gap+nSFA), comp2, zeros(2,ngtime_gap), ...
    SFA(6)*PSFA0, zeros(2,ngtime_gap), comp1, zeros(2,2*ngtime_gap+nSFA), comp1, zeros(2,ngtime_gap), ...
    SFA(7)*PSFA0, zeros(2,ngtime_gap), comp1, zeros(2,2*ngtime_gap+nSFA), comp2, zeros(2,ngtime_gap), ...
    SFA(8)*PSFA0, zeros(2,ngtime_gap), comp2, zeros(2,2*ngtime_gap+nSFA), comp1, zeros(2,ngtime_gap), ...
    SFA(9)*PSFA0 ];

% add additional crusher gradient
g_crusher = Ggen('t',[1],0.32*[1],3*[1],[0.08, 0.28],dt); % 0.08 and 0.28 are gaps (ms) before and after the crusher
ng_crusher = numel(g_crusher);
g = [g, g_crusher];
rf = [rf, zeros(2,ng_crusher)];

% limit phase to be within -pi to pi
rf(2,:) = mod(rf(2,:)+pi,2*pi)-pi;

%% plot pulses
figure;
plot(rf(1,:),'r');
hold on;
plot(rf(2,:),'g');
plot(g,'b');
hold off;

% %% output the waveform as txt file
% % rho
% nrf = size(rf,2);
% fileID = fopen(['myVSI_' num2str(nrf) '.rho.txt'],'w'); 
% for i=1:nrf
%     num = 2*round(rf(1,i)*16383);
%     fprintf(fileID,'%i\n',num);
% end
% fclose(fileID);
% 
% % theta
% fileID = fopen(['myVSI_' num2str(nrf) '.theta.txt','w'); 
% for i=1:nrf
%     num = 2*round(rf(2,i)/pi*16383);
%     fprintf(fileID,'%i\n',num);
% end
% fclose(fileID);
% 
% % grad
% fileID = fopen(['myVSI_' num2str(nrf) '.grad.txt','w'); 
% for i=1:nrf
%     num = 2*round(g(i)*16383);
%     fprintf(fileID,'%i\n',num);
% end
% fclose(fileID);
