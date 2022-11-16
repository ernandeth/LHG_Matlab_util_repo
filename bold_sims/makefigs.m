


%%%%%%%%%%%%%%%%
% plots of HRF errors
load tau_error
sim = avg_t;

figure
t = [-1:0.1:1];
plot(t,avg_t(1,:),'--k','linewidth',3)
hold on
plot(t,avg_t(11,:),'k','linewidth',3)
plot(t,avg_t(21,:),'o-k','linewidth',3)

title('HRF errors', 'fontsize',16, 'fontweight','bold')
xlabel('Temporal Shift of the Reggressors (seconds)', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('No Error','Error in Tau=0.5', 'Error in Tau=-0.5',4)

print -dtiff newfigs/tau_error.tif


load tau_error_noisy
sim = avg_t;

figure
t = [-1:0.1:1];
plot(t,sim(1,:),'--k','linewidth',3)
hold on
plot(t,sim(11,:),'k','linewidth',3)
plot(t,sim(21,:),'o-k','linewidth',3)

title('HRF errors (noise)', 'fontsize',16, 'fontweight','bold')
xlabel('Temporal Shift of the Reggressors (seconds)', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('No Error','Error in Tau=0.5', 'Error in Tau=-0.5',4)

print -dtiff newfigs/tau_error_noisy.tif


load tau_error_data
dat = avg_t;

figure
t = [-1:0.1:1];
plot(t,dat(1,:),'--k','linewidth',3)
hold on
plot(t,dat(11,:),'k','linewidth',3)
plot(t,dat(21,:),'o-k','linewidth',3)
title('HRF errors (data)', 'fontsize',16, 'fontweight','bold')
xlabel('Temporal Shift of the Regressor (seconds)', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('Error in Tau = -0.5', 'No Error', 'Error in Tau = 0.5',4)

print -dtiff newfigs/tau_error_data.tif



figure

dt = [-0.5:0.05:0.5];
plot(dt, dat(:,11),'k','linewidth',3)
hold on 
plot(dt, sim(:,11), '--k','linewidth',3);


title('HRF errors (no shift)', 'fontsize',16, 'fontweight','bold')
xlabel('Error in the Tau constants (seconds)', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('Data','Simulation',4)

print -dtiff newfigs/tau_error_noshift.tif



%%%%%

% plots of the ITI errors:

load iti_error
sim = avg_t;

figure
t = [-1:0.1:1];
plot(t,sim(1,:),'k','linewidth',3)
hold on
plot(t,sim(11,:),'--k','linewidth',3)

title('ITI errors', 'fontsize',16, 'fontweight','bold')
xlabel('Temporal Shift of the Regressor (seconds)', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('No ITI Error', 'ITI error jitter = 1',4)

print -dtiff newfigs/iti_error.tif

load iti_error_noisy
sim = avg_t;

figure
t = [-1:0.1:1];
plot(t,sim(1,:),'k','linewidth',3)
hold on
plot(t,sim(11,:),'--k','linewidth',3)

title('ITI errors (noisy)', 'fontsize',16, 'fontweight','bold')
xlabel('Temporal Shift of the Regressor (seconds)', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('No ITI Error', 'ITI error jitter = 1',4)

print -dtiff newfigs/iti_error_noisy.tif


load iti_error_data
dat = avg_t;
figure

t = [-1:0.1:1];
plot(t,dat(1,:),'k','linewidth',3)
hold on
plot(t,dat(11,:),'--k','linewidth',3)

title('ITI errors (data)', 'fontsize',16, 'fontweight','bold')
xlabel('Temporal Shift of the Regressor (seconds)', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('No ITI Error', 'ITI error jitter = 1',4)

print -dtiff newfigs/iti_error_data.tif


figure
d_iti = [0:0.1:1];
plot(d_iti, dat(:,11),'k','linewidth',3)
hold on 
plot(d_iti, sim(:,11), '--k','linewidth',3);

title('ITI errors (no shift)', 'fontsize',16, 'fontweight','bold')
xlabel('Variance of the ITI Error (seconds)', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('Data','Simulation',4)

print -dtiff newfigs/iti_error_noshift.tif

%%%%%%%%%%%%%
load n_iti_error
dat = avg_t;

figure
t = [-1:0.1:1];
plot(t,avg_t(1,:),'o-k','linewidth',3)
hold on
plot(t,avg_t(11,:),'k','linewidth',3)
plot(t,avg_t(21,:),'--k','linewidth',3)
title('Wrong number of trials', 'fontsize',16, 'fontweight','bold')
xlabel('Temporal Sensitivity (seconds)', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('Number of Extra Events = -10', ...
   'perfect match', ...
   'Number of extra events = 10',4)

print -dtiff newfigs/niti_error.tif


load n_iti_error_data
dat = avg_t;
load n_iti_error_noisy
sim = avg_t;


figure
t = [-1:0.1:1];
plot(t,avg_t(1,:),'o-k','linewidth',3)
hold on
plot(t,sim(11,:),'k','linewidth',3)
plot(t,sim(21,:),'--k','linewidth',3)
title('Wrong number of trials  (noisy)', 'fontsize',16, 'fontweight','bold')
xlabel('Temporal Sensitivity (seconds)', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('Number of Extra Events = -10', ...
   'perfect match', ...
   'Number of extra events = 10',4)

print -dtiff newfigs/niti_error_noisy.tif

figure
t = [-1:0.1:1];
plot(t,dat(1,:),'o-k','linewidth',3)
hold on
plot(t,dat(11,:),'k','linewidth',3)
plot(t,dat(21,:),'--k','linewidth',3)
title('Wrong number of trials (data)', 'fontsize',16, 'fontweight','bold')
xlabel('Temporal Sensitivity (seconds)', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('Number of Extra Events = -10', ...
   'perfect match', ...
   'Number of extra events = 10',4)

print -dtiff newfigs/niti_error_data.tif

figure
d_niti = [-10:10];
plot(d_niti, dat(:,11),'k','linewidth',3)
hold on 
plot(d_niti, sim(:,11), '--k','linewidth',3);
ax = axis;
ax(1) = -10;
axis(ax);

title('Wrong number of trials (no shift)', 'fontsize',16, 'fontweight','bold')
xlabel('Number of Errors', 'fontsize',16, 'fontweight','bold')
ylabel('t score', 'fontsize',16, 'fontweight','bold')
set(gca,'fontsize',12, 'fontweight','bold','linewidth',2)
legend('Data','Simulation',4)

print -dtiff newfigs/niti_error_noshift.tif

return
%%%%%%%%%%%%%%%%
load spm_vs_myglm.mat

plot(t_spm(:,1), t_spm(:,2),'--k','linewidth',3);
hold on
plot(t(:,1), t(:,2),'k','linewidth',3);

title('Visual Stimulation Experiment Time-shift Analysis')
xlabel('shift')
ylabel('t score')

legend('single pix (SPM)', '5mm Sphere (my-glm)',4)
