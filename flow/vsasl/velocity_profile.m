% velocity selective profile simulator
% for BIR-8 pulse , we use equation 2 in Meakin 2013
% this is the average Mz returned to the z axis for a parabolic velocity profile


gamma = 2*pi*4.257e3; % rad/s/gauss
% First moment of the gradient in the BIR8 pulses 
% (2ms triangle with 1 ms gap after it)
% run the BI_sim program to construct the pulse sequence
vel=20
bir_sim_06

%%
% Using the equation from the Meakin paper (2013)
gamma = 2*pi*4.257e3; % rad/s/gauss
sampTime=4e-6;
t=[0:length(Gz)-1] * sampTime;
Vmax_range = linspace(0,10,500);

m1 = sum(t .* abs(Gz))*sampTime;

gz=load('pulses/myVSI_6850.grad.txt');
gz = 1*gz(1:6000)/max(gz);
gz = gz';
% to get a cutoff near 1 cm/s you need about 2 G/cm Gmax
m1 = sum(t(1:6000) .* abs(gz))*sampTime;

profile = sinc(gamma*m1*Vmax_range);

% cutoff velocity is when Mz crosses zero
c_ind = find(profile<0);
Vcutoff_BIR = Vmax_range(c_ind(1))
fprintf('\nOfficial Vcut (first zero crossing) for BIR8 is %f', Vcutoff_BIR);

figure
plot(Vmax_range, profile), grid on
axis([0 max(Vmax_range) -1 1])
xlabel('Velocity (cm/s)')
ylabel('M_z (a.u.)')
title('BIR8 Velocity  Profile')
fatlines
dofontsize(16)
legend off
print -dpng ~/BIR8_profile

%%  now do the same thing by brute numerical force
%{
all_Mz=[];
Vmax_range = linspace(0,80,150);

for vel=Vmax_range
    bir_sim_06
    Mz = M(end,3);
   all_Mz = [all_Mz Mz];
end
dv = Vmax_range(2)-Vmax_range(1);



% cutoff velocity
c_ind = find(profile<0.1);
Vcutoff = Vmax_range(c_ind(1)-1)
%}

%% We will now do the paraboilc profile with the FTVS pulses
%vel_sim_09
%%
% vel_sim_11  % for the 17536 pulse - FTVSS
vel_sim_12  
figure
plot(vel_range*1e3, Mzfinal)
hold on
%plot(vel_range*1e3, Mzfinal_ns)
axis([0 10 0 1])
grid on
%legend('Selective', 'Non-Selective')
xlabel('Velocity (cm/s)')
ylabel('M_z (a.u.)')
title('FTVS Profile')
fatlines
dofontsize(16)
legend off
print -dpng ~/FTVS_profile


% the MEAN velocity profile should be determined by the sum of the Mz 
% from all the velocities in the range of velocities 
% from 0 to Vmax.
% weighted by the profile (parabola peaking at Vpeak).
%
%delta = ones(size(delta));
profile = zeros(size(delta));
for v=1:length(vel_range)
    Vp = vel_range(v);
    weight = [0:length(vel_range)-1]';
    weight = -(weight-v).^2 + (v)^2;
    wm = find(weight == max(weight));
    weight(wm+1:end) = 0;
    weight = weight/sum(weight);
    
    profile(v) = sum(delta(1:v) .* weight(1:v));
    
end
figure
plot(vel_range*1e3, profile)
title('Mean velocity profile for FTVS pulse')
xlabel('Velocity (cm/s)')
ylabel('M_z (a.u.)')


%%
t2 = [0:length(Gz)-1]'*1e-6;
Gz = Gz/max(Gz);
m1_ftvs = sum(t2 .* abs(Gz));
gamma = 2*pi*4.257e3; % rad/s/gauss
profile_ftvs = sinc(gamma*m1_ftvs*Vmax_range);
figure
plot(Vmax_range, profile_ftvs), grid on



