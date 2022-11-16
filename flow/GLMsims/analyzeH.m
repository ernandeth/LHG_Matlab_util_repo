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

betas = [1e4  1e2 1.5*1e2 1e2]';

for ypix=1:h.ydim
    for xpix=1:h.xdim

        noise = (NoiseLevel) * randn(size(Nframes,1));
        ind = sub2ind([h.xdim, h.ydim], xpix, ypix);
        Y(:,ind) = Dfull*betas + noise;
        %plot(Y(ind,:)); drawnow
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

beta2flow02('ConBhats','ConVar_hats', TR,Ttag,pid,Ttrans, inv_alpha, 0);

flows = read_img('ExpFlows');

f0 = reshape(flows(2,:), h.xdim, h.ydim, 1); imagesc(f0); colorbar
f1 = reshape(flows(3,:), h.xdim, h.ydim, 1) ;imagesc(f1); colorbar

fvars = read_img('ExpFlow_vars');
fv0 = reshape(fvars(2,:), h.xdim, h.ydim, 1); imagesc(fv0); colorbar
fv1 = reshape(fvars(3,:), h.xdim, h.ydim, 1) ;imagesc(fv1); colorbar


%%
%%% Analysis 2 - Traditional Way

raw = Y;

rest=[	4:11		31:38		56:63		81:88] + 1;
active = [	18:25		43:50		68:75		93:100];

rest_tag = raw(rest(1:2:end),:);
rest_con = raw(rest(2:2:end),:);

act_tag = raw(active(1:2:end),:);
act_con = raw(active(2:2:end),:);

h.tdim = size(rest_tag,1);

write_img('rest_tag.img',abs(rest_tag),h);
write_img('rest_con.img',abs(rest_con),h);
write_img('act_tag.img',abs(act_tag),h);
write_img('act_con.img',abs(act_con),h);

[act_ff  act_ffv] = calcPerf_v('act_con', 'act_tag', 'actPerf.img', TR,Ttag,pid,Ttrans, inv_alpha);
[rest_ff rest_ffv] = calcPerf_v('rest_con', 'rest_tag', 'restPerf.img', TR,Ttag,pid,Ttrans, inv_alpha);


f0 = reshape(rest_ff, h.xdim, h.ydim, 1); imagesc(f0); colorbar
f1 = reshape(act_ff, h.xdim, h.ydim, 1) ;imagesc(f1); colorbar

fv0 = reshape(rest_ffv, h.xdim, h.ydim, 1); imagesc(fv0); colorbar
fv1 = reshape(act_ffv, h.xdim, h.ydim, 1) ;imagesc(fv1); colorbar

%% show what was averaged

plot(Dfull(:,2:3))
hold on; plot(rest, Dfull(rest,2), '*');

hold on; plot(active, Dfull(active,3), 'g*');

axis([0 30 -1.2 0.8])
legend('Resting ASL signal', 'Active ASL signal', 'Resting Samples', 'Active Samples')
title('Sampling for Perfusion Calculation')
xlabel('scan number')
ylabel('A.U.')
fatlines
dofontsize(12)