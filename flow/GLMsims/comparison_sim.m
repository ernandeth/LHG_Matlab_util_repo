%% simulation code for testing the perfusion estimator using general linear
%% model parameter estimates
%% First Build the Design Matrix
TR = 4;
Ttag = 2.2;
pid = 1.5;
Ttrans =1.5 ;
inv_alpha = 0.8;

lambda = 0.9;
T1 = 1.25;
T1a = 1.65;
T1app=1/(1/T1 + 0.015/lambda);

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
    %D( rest_on(c) + 2*TR: (rest_on(c) + rest_dur(c) ) , 3) = 1;
    D(: , 3) = 1;
end



h = spm_hrf(1);
for r=2:2
    reg = D(:,r);
    reg = conv(reg,h);
    reg = reg(1:duration);
    D(:,r) =  reg;
end

% downsample:
D = D(1:TR:end, :);

D2 = D;

% ASL modulation of regressors
D2(1:2:end,2:end) = -D2(1:2:end,2:end);
D2 = D2/2;

Dfull= zeros(abs(duration/TR) , 4);
Dfull (:,1) = D(:,1);
Dfull(:,3) = D2(:,2);
Dfull(:,2) = D2(:,3);
Dfull(:,4) = D(:,2);%-mean(D(:,2));

figure
imagesc(Dfull)
colormap(gray)
xlabel('Regressor')
ylabel('Scan Number')
title('Design Matrix')
set(gca,'XTick',[ 1 2 3 4]);
%'\beta_0' , '\beta_1', '\beta_2', '\beta_3'])
dofontsize(16)

figure
subplot(411)
plot(Dfull(:,1), 'k'), title('Baseline MR Signal (x _{0t})')
dofontsize(16); 
subplot(412)
plot(Dfull(:,2), 'k'), title('Baseline ASL Signal (x _{1t})')
 dofontsize(16); axis tight; 
subplot(413)
plot(Dfull(:,3), 'k'), title('Activation ASL Signal (x _{2t})')
dofontsize(16); axis tight ; 
subplot(414)
plot(Dfull(:,4), 'k'), title('Activation BOLD Signal (x _{3t})')
xlabel('Scan Number')
dofontsize(16); axis tight; 
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

betas = [1e4  0.5e2 0.2e2 0.5e2]';
trubetas = betas;

% I'm going to take advantage of existing code that works on images...  
% all pixels along the x dimension are instances of the same noise level
% the noise increases along the y dimension
% there is only one slice

for ypix=1:h.ydim
    for xpix=1:h.xdim

        noise = sqrt(NoiseLevel) * randn(Nframes,1);
        ind = sub2ind([h.xdim, h.ydim], xpix, ypix);
        Y(:,ind) = Dfull*betas + noise;
        %plot(Y(:, ind)); drawnow
    end
    NoiseLevel = NoiseLevel + 10;
end
write_img('fakedata.img', Y, h);


%% calculation of "true perfusions" given the true betas (this is what's in
% the main program
Nbetas = size(betas,1);

M0 = betas(1)/(1 - exp(-TR/T1));

% This is equation 3 from Alsop et al: JCBFM 16, 1236-1249,1996
%den = T1app*2*M0*inv_alpha/lambda * exp(-Ttrans*(1/T1a-1/T1app))*exp(-pid/T1a);

% modification by Wang et al MRM 48,2,p242-254, 2002:
%den = T1app*2*M0*inv_alpha/lambda * ...
%    (exp(-Ttrans*(1/T1a-1/T1app)) - exp((Ttrans-pid-Ttag)/T1app))...
%    *exp(- pid/T1a);

% looking at equation 1, delta_a and delta should be almost the same, we
% can drop a term.
den =  2 * M0* (inv_alpha / lambda)...
    * T1app * exp(-Ttrans/T1a)...
	* ( exp( (Ttrans-pid)/T1app ) - exp( (Ttrans-Ttag-pid)/T1app));

den = repmat(den, Nbetas,1);
f = betas./den * 6000;  % conversion factor from ml/s/g to ml/min/100g






%%  Analysis 1 _ GLM and estimation
flags.header = h;
flags.doWhiten = 0;

spmJr(Y , Dfull, ...
    [...
    1 0 0 0;...
    0 1 0 0;
    0 0 1 0;
    0 0 0 1],flags);

close all

beta2flow02('ConBhats','ConVar_hats', TR, Ttag, pid, Ttrans, inv_alpha, 0);

flows = read_img('ExpFlows');

figure

f0 = reshape(flows(2,:), h.xdim, h.ydim, 1); imagesc(f0); colorbar
f1 = f0+reshape(flows(3,:), h.xdim, h.ydim, 1) ;imagesc(f1); colorbar

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


tf0 = -reshape(rest_ff, h.xdim, h.ydim, 1); imagesc(tf0); colorbar
tf1 = reshape(act_ff, h.xdim, h.ydim, 1) ;imagesc(tf1); colorbar

tfv0 = reshape(rest_ffv, h.xdim, h.ydim, 1); imagesc(tfv0); colorbar
tfv1 = reshape(act_ffv, h.xdim, h.ydim, 1) ;imagesc(tfv1); colorbar

%% show what was averaged
figure
subplot(211)


hold on; plot(rest, Dfull(rest,2), 'kx');
title('Sampling for Perfusion Calculation')
xlabel('Scan Number')
ylabel('Rest')
fatlines
plot(Dfull(:,2),'k');
axis tight
dofontsize(12)

subplot(212)
hold on
plot(active, Dfull(active,3), 'kx');
fatlines
plot(Dfull(:,3),'k')
axis tight
ylabel('Active')
xlabel('Scan Number')

dofontsize(12)


%%  Show the correlation between the measuremens
figure
subplot(121), plot(mean(f1,1),mean(tf1,1),'*'), axis square, 
axis([54 56 54 56])
title('Active '), xlabel('GLM'), ylabel('Traditional');
r =corrcoef(mean(f1,1), mean(tf1,1))
text( 54, 55, sprintf('R = %0.2f', r(1,2) ))
fatlines,dofontsize(12)

subplot(122), plot(mean(f0,1),mean(tf0,1),'o'), axis square, 
axis([38 42 38 42])

title('Rest'), xlabel('GLM'), ylabel('Traditional');
r2 =corrcoef(mean(f0,1), mean(tf0,1))
text(39,40, sprintf('R = %0.2f', r2(1,2) ))
fatlines,dofontsize(12)

x = mean(tf0,1)';
y = mean(f0,1)';
Xmat=[zeros(size(y)) x];
betas0 = pinv(Xmat)*y

x = mean(tf1,1)';
y = mean(f1,1)';
Xmat=[zeros(size(y)) x];
betas1 = pinv(Xmat)*y

%% show the different variance estimates
figure
subplot(121), 
plot(sqrt([0:10:499]), sqrt(mean(tfv0,1)), 'k')
hold on, 
plot(sqrt([0:10:499]),sqrt(mean(tfv1,1)),'--k')
legend('Rest', 'Active')
title('Traditional')
fatlines; dofontsize(12)
xlabel('Added Noise Std. Dev. (a.u.)');
ylabel(' Perf. Std. Dev. (ml/min/100g)');
fatlines; dofontsize(12)

subplot(122), 
plot(sqrt([0:10:499]), sqrt(mean( fv0,1)), 'k')
hold on, 
plot(sqrt([0:10:499]), sqrt(mean(fv1,1)),'--k')
legend('Rest', 'Active')
title('GLM')

xlabel('Added Noise Std. Dev. (a.u.)');
ylabel(' Perf. Std. Dev. (ml/min/100g)');
fatlines; dofontsize(12)


%% Now let's see how the experimental data did: human data
Res = load('Results26-Aug-2009');

r = Res.Results

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
xlabel('GLM estimates (ml/min/100g)')
ylabel('Traditional Model Means (ml/min/100g)')
axis square; axis([30 90 30 90])
rho = corrcoef(f,ft)
text(40,80, sprintf('R = %0.2f', rho(1,2)));
fatlines; dofontsize(12)
 
 
%% timeseries plots:  pairwise subtraction vs. estimation

% first plot the truth
Nbetas = size(trubetas,1);
M0 = trubetas(1)/(1 - exp(-TR/T1));
den =  2 * M0* (inv_alpha / lambda)...
    * T1app * exp(-Ttrans/T1a)...
	* ( exp( (Ttrans-pid)/T1app ) - exp( (Ttrans-Ttag-pid)/T1app));

den = repmat(den, Nbetas,1);
f = trubetas./den * 6000;  % conversion factor from ml/s/g to ml/min/100g

% the relevant fs are the 2nd and third regressors
truflow = 2*Dfull(:,2:3) * f(2:3);
truflow(1:2:end) = -truflow(1:2:end);

figure
plot(truflow, 'k');
hold on

mypix = 500; % (corresponds to NoiseLevel = (mypx-1)*10)
% pairwise subtraction
data = Y(:, mypix);
[m v f] = calcPerf_v(data(2:2:end), data(1:2:end), 'restPerf.img', TR,Ttag,pid,Ttrans, inv_alpha);

t = [1:length(f)]*2;
plot(t,f,'g');
title('Perfusion timeseries')

rest_f = f(floor(rest/2));
act_f = f(floor(active/2));

% traditional results
tradf_0 = mean(rest_f);
tradf_1 = mean(act_f);
tradfv_0 = var(rest_f);
tradfv_1 = var(act_f);

% GLM:
allestflows = read_img('ExpFlows.img');
allbetas = read_img('conBhats.img');
allbetavars = read_img('ConVar_hats');
allvars = read_img('expFlow_vars');

% GLM results
betas = allbetas(:,mypix);
vars = allvars(:,mypix);
betavars = allbetavars(:,mypix);
f = allestflows(:,mypix);

% Nbetas = size(betas,1);
% M0 = betas(1)/(1 - exp(-TR/T1));
% den =  2 * M0* (inv_alpha / lambda)...
%     * T1app * exp(-Ttrans/T1a)...
% 	* ( exp( (Ttrans-pid)/T1app ) - exp( (Ttrans-Ttag-pid)/T1app));
% 
% den = repmat(den, Nbetas,1);
% f = betas./den * 6000;  % conversion factor from ml/s/g to ml/min/100g

% the relevant fs are the 2nd and third regressors
estflow = 2*Dfull(:,2:3) * f(2:3);
estflow(1:2:end) = -estflow(1:2:end);


%estflow = estflow(1:2:end) - estflow(2:2:end);
hold on
plot(estflow )

xlabel('Scan Pair')
ylabel('Perfusion (ml/min/100g)')
axis ([ 1 100 0 100])
legend('True Perfusion','Traditional Method', 'GLM estimated')
dofontsize(10)
fatlines
plot(estflow -sqrt(vars(2)+vars(3)),'--')
plot(estflow +sqrt(vars(2)+vars(3)),'--')



%% get some time courses from experimental human data:
data = load('vis_rawtdata.dat')
[m v f] = calcPerf_v(data(2:2:end), data(1:2:end), 'restPerf.img', TR,Ttag,pid,Ttrans, inv_alpha);

t = [1:length(f)]*2;
figure
plot(t,f,'g');
title('Perfusion timeseries (Human Data)')


% traditional results
rest_f = f(floor(rest/2));
act_f = f(floor(active/2));

tradf_0 = mean(rest_f);
tradf_1 = mean(act_f);
tradfv_0 = var(rest_f);
tradfv_1 = var(act_f);

% GLM:
flags.header = [];
flags.doWhiten = 1;



[betas bvars ] =spmJr(data , Dfull, ...
    [...
    1 0 0 0;...
    0 1 0 0;
    0 0 1 0;
    0 -1 1 0],flags);


Nbetas = size(betas,1);
M0 = betas(1)/(1 - exp(-TR/T1));
den =  2 * M0* (inv_alpha / lambda)...
    * T1app * exp(-Ttrans/T1a)...
	* ( exp( (Ttrans-pid)/T1app ) - exp( (Ttrans-Ttag-pid)/T1app));

den = repmat(den, Nbetas,1);
f = betas./den * 6000;  % conversion factor from ml/s/g to ml/min/100g

% the relevant fs are the 2nd and third regressors
estflow = 2*Dfull(:,2:3) * f(2:3);
estflow(1:2:end) = -estflow(1:2:end);


% now the vars

varM0 = bvars(1);

df_dBeta = lambda * (1/T1app) ./ ...
	(...
	M0 * 2 * inv_alpha ...
    * T1app * exp(-Ttrans/T1a) ...
	* ( exp( (Ttrans-pid)/T1app ) - exp( (Ttrans-Ttag-pid)/T1app))...
	);

df_dM0 = - betas * lambda * (1/T1app) ./ ...
	(...
	(M0.^2) * 2 * inv_alpha ...
    * T1app * exp(-Ttrans/T1a) ...
	* ( exp( (Ttrans-pid)/T1app ) - exp( (Ttrans-Ttag-pid)/T1app))...
	);

flowVars =   varM0 .* df_dM0.^2  + vars .* df_dBeta.^2;


flowVars = flowVars * 6000^2;  % units conversion

%estflow = estflow(1:2:end) - estflow(2:2:end);
hold on
plot(estflow )

xlabel('Scan Pair')
ylabel('Perfusion (ml/min/100g)')
axis ([ 1 100 0 140])
legend('Traditional Method', 'GLM estimated')
dofontsize(10)
fatlines
plot(estflow -sqrt(flowVars(2)+flowVars(3)),'--')
plot(estflow +sqrt(flowVars(2)+flowVars(3)),'--')




