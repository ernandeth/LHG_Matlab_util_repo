%% simulation code for testing the perfusion estimator using general linear
%% model parameter estimates
%% First Build the Design Matrix
TR = 4;
Ttag = 2.2;
pid = 1.5;
Ttrans =1.5 ;
inv_alpha = 0.8;

Nframes = 100;
duration = Nframes*TR;


%%%%% Activation
act_on = [50:100:duration];
rest_on = [1:100:duration];

%%%%% Rest
act_dur = [50 50 50 50];
rest_dur = [50 50 50 50];


D = zeros(duration , 3);
%\beta0

%\beta 0
D(:,1)=1;

%\beta1:rest beta2:act perfusion beta3: act BOLD
for c=1:length(act_on)
    D( act_on(c)  : (act_on(c) + act_dur(c) ) , 2) = 1;
    D( rest_on(c) : (rest_on(c) + rest_dur(c) ) , 3) = 1;
end



h = spm_hrf(1);
for r=2:3
    reg = D(:,r);
    reg = conv(reg,h);
    reg = reg(1:duration);
    D(:,r) =  reg;
end

D = D(1:TR:end, :);
D2 = D;
D2(1:2:end,:) = -D2(1:2:end,:);


Dsub = D(1:2:end,:);
% surround subtraction case
%Dsursub = D(2:end-1, :);
Dsursub = -differencer(D2,4);

D2 = D2/2;
imagesc([D D2])

Dfull= zeros(abs(duration/TR) , 4);
Dfull (:,1)=D(:,1);
Dfull(:,3)=D2(:,2);
Dfull(:,2)=D2(:,3);
Dfull(:,4)=D(:,2);
%%

% generate data using the above matrix
h = define_avw_hdr;
h.xdim = 100;
h.ydim = 50;
h.zdim = 1;
h.tdim = 100;
h.dims = 4;
h.xsize = 1;
h.ysize = 1;
h.zsize = 1;
h.datatype = 8;
h.bits = 16;

Y = zeros(h.tdim, h.xdim*h.ydim);
NoiseLevel=0;

betas = [1e4  1e2 1.2*1e2 1e2]';

for ypix=1:h.ydim
    for xpix=1:h.xdim

        noise = sqrt(NoiseLevel) * randn(Nframes,1);
        ind = sub2ind([h.xdim, h.ydim], xpix, ypix);
        Y(:,ind) = Dfull*betas + noise;
        %plot(Y(:, ind)); drawnow
    end
    NoiseLevel = NoiseLevel + 1;
end


write_img('fakedata.img', Y, h);

%%  Analysis 1 _ GLM and estimation
flags.header = h;
flags.doWhiten = 0;

spmJr(Y , Dfull, ...
    [...
    1 0 0 0;...
    0 1 0 0;
    0 0 1 0;
    0 -1 1 0],flags);

close all

beta2flow02('ConBhats','ConVar_hats', TR, Ttag, pid, Ttrans, inv_alpha, 0);

flows = read_img('ExpFlows');

f0 = reshape(flows(2,:), h.xdim, h.ydim, 1); imagesc(f0); colorbar
f1 = reshape(flows(3,:), h.xdim, h.ydim, 1) ;imagesc(f1); colorbar

fvars = read_img('ExpFlow_vars');
fv0 = reshape(fvars(2,:), h.xdim, h.ydim, 1); imagesc(fv0); colorbar
fv1 = reshape(fvars(3,:), h.xdim, h.ydim, 1) ;imagesc(fv1); colorbar


%%
%%% Analysis 2 - Traditional Way

raw = Y;

rest=[	5:12		31:38		57:64		81:88] ;
active = [	19:26		43:50		69:76		93:100];

rest_con = raw(rest(1:2:end),:);
rest_tag = raw(rest(2:2:end),:);

act_tag = raw(active(1:2:end),:);
act_con = raw(active(2:2:end),:);

h.tdim = size(rest_tag,1);

write_img('rest_tag.img',abs(rest_tag),h);
write_img('rest_con.img',abs(rest_con),h);
write_img('act_tag.img',abs(act_tag),h);
write_img('act_con.img',abs(act_con),h);

[act_ff  act_ffv] = calcPerf_v('act_con', 'act_tag', 'actPerf.img', TR,Ttag,pid,Ttrans, inv_alpha);
[rest_ff rest_ffv] = calcPerf_v('rest_con', 'rest_tag', 'restPerf.img', TR,Ttag,pid,Ttrans, inv_alpha);


tf0 = reshape(rest_ff, h.xdim, h.ydim, 1); imagesc(tf0); colorbar
tf1 = reshape(act_ff, h.xdim, h.ydim, 1) ;imagesc(tf1); colorbar

tfv0 = reshape(rest_ffv, h.xdim, h.ydim, 1); imagesc(tfv0); colorbar
tfv1 = reshape(act_ffv, h.xdim, h.ydim, 1) ;imagesc(tfv1); colorbar

%% show what was averaged
figure
subplot(211)

axis([0 30 -1.2 0.8])
hold on; plot(rest, Dfull(rest,2), 'kx');
title('Sampling for Perfusion Calculation')
xlabel('scan number')
ylabel('Rest')
fatlines
plot(Dfull(:,2),'k')
dofontsize(12)

subplot(212)
hold on
plot(active, Dfull(active,3), 'kx');
fatlines
plot(Dfull(:,3),'k')
axis([0 30 -1.2 0.8])
ylabel('Active')
xlabel('scan number')

dofontsize(12)


%%  Show the correlation between the measuremens
figure
subplot(121), plot(mean(f1,1),mean(tf1,1),'*'), axis square, axis([94 95 92 93])
title('Active '), xlabel('GLM'), ylabel('Traditional');
r =corrcoef(mean(f1,1), mean(tf1,1))
text( 94.1, 92.8, sprintf('R = %0.2f', r(1,2) ))
fatlines,dofontsize(12)

subplot(122), plot(mean(f0,1),mean(tf0,1),'o'), axis square, axis([78 79 75.5  76.5])
title('Rest'), xlabel('GLM'), ylabel('Traditional');
r =corrcoef(mean(f0,1), mean(tf0,1))
text(78.1,76.2, sprintf('R = %0.2f', r(1,2) ))
fatlines,dofontsize(12)


%% show the different variance estimates
figure
subplot(121), 
plot([0:49], mean(tfv0,1))
hold on, 
plot([0:49],mean(tfv1,1),'r')
legend('Rest', 'Active')
title('Traditional')
fatlines; dofontsize(12)
xlabel('Noise Variance (a.u.)');
ylabel('Variance of Perf. Measurement (ml/min/100g');

subplot(122), 
plot([0:49], mean( fv0,1))
hold on, 
plot([0:49], mean(fv1,1),'r')
legend('Rest', 'Active')
title('GLM')
fatlines; dofontsize(12)
xlabel('Noise Variance (a.u.)');
ylabel('Variance of Perf. Parm. Estiamtes (ml/min/100g');

%% Now let's see how the experimental data did:
Res = load('RESULTS_SMALL_ROIS');
r = Res.Results;

% concatenate all the motor active and rest
% GLM estimates
f = [r(:,1); % motor act 
    r(:,3);  % motor base
    r(:,5);  % visual act
    r(:,7);  % visual base];
    ];
% traditional
ft = [r(:,9); % motor act 
    r(:,11);  % motor base
    r(:,13);  % visual act
    r(:,15);  % visual base];
    ];

figure
plot(f,ft,'o')
title('Comparison of Perfusion Estimates')
xlabel('GLM estimatess (ml/min/100g')
ylabel('Traditional Model Means (ml/min/100g)')
axis square; axis([30 90 30 90])
fatlines; dofontsize(12)

