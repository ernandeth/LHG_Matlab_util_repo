

%%
Nvels = 50;
all_maxVels = linspace(0,30,Nvels);

Mz_batch = zeros(Nvels,1);
M_final = zeros(Nvels, 3);
ph_final = zeros(Nvels,1);

m=1;

for maxVel= all_maxVels
    vels = linspace(0, maxVel,Nvels);
    
    v = 1;
    for vel=vels
        % LHG 3/23/22
        % bir_sim_06
        %
        bir_sim_lofi_220323
        M_final(v,:)    = M(end,:);
        v = v+1;
    end
    
    signal(m) = sum(M_final(:,3))/Nvels;
    m = m+1
end

eTE = -log(M_final(1,3))*T2
%{
figure(1);
subplot (311)
plot(vels,M_final(:,1)), hold on
plot(vels,M_final(:,2))
xlabel('velocity')
ylabel('M_{xy}')
legend('M_x', 'M_y')
grid on

subplot (312)
plot(vels,M_final(:,3)), hold on
xlabel('velocity')
ylabel('M_z')

subplot(313)
plot(vels, signal)
xlabel('velocity')
ylabel('Net M_z')
grid on
%}

%%
figure(2)
hold on
plot(vels, signal)
xlabel('Velocity (cm/s)')
ylabel('Net  M_z')
grid on
title('BIR-8 VSS Velocity Profiles')