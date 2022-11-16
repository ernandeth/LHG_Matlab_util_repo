function calc_mrf_errors_20220225

vfiles = dir('vol*.nii');
[raw h]= read_nii_img([vfiles(1).name]);

[cbf h] = read_nii_img('CBF_test.nii');
cbv = read_nii_img('CBV_test.nii');
bat = read_nii_img('BAT_test.nii');
t1 = read_nii_img('T1_test.nii');
t2 = read_nii_img('T2_test.nii');
flip = read_nii_img('flip_test.nii');

bat(bat<0) = 0.01;

error = zeros(size(cbf));
Npix = length(error(:));


formatSpec = '%f';
fileID = fopen('t_delays.txt','r');
timing_parms.t_delay = fscanf(fileID,formatSpec);
fileID = fopen('t_adjusts.txt','r');
timing_parms.t_adjusts = fscanf(fileID,formatSpec);
fileID = fopen('doArtSuppression.txt','r');
timing_parms.doArtSup =  fscanf(fileID,formatSpec);
fileID = fopen('isVelocitySelective.txt','r');
timing_parms.labelcontrol =  fscanf(fileID,formatSpec);

% some preset constant values:
Nframes = length(timing_parms.t_delay);
timing_parms.t_tag  = 0.1*ones(Nframes,1) ;
timing_parms.ArtSup_delay  = 0.15*ones(Nframes,1) ;
timing_parms.t_aq  = 0.0329*18*ones(Nframes,1) ;
timing_parms.label_type = 'FTVSI-sinc';
timing_parms.readout_type = 'FSE';
timing_parms.order = 1;

%
parms=[];
parms.f     = 0;
parms.cbva  = 0;
parms.bat   = 0;
parms.r1tis = 0;
parms.r2tis = 0;
parms.flip  = deg2rad(45);

parms.mtis0 =     1 ;
parms.Disp =      40;
parms.r1blood =   1/1.7;

%parms = repmat(parms,Npix,1);
parms = repmat(parms, Npix,1);

Nframes = 200;
raw = raw(1:Nframes, :);
cooked = zeros(size(raw));

parfor p=1:Npix
    parms(p).f     = cbf(p)/6000;
    parms(p).cbva  = cbv(p);
    parms(p).bat   = bat(p);
    parms(p).r1tis = 1/t1(p);
    parms(p).r2tis = 1/t2(p);
    parms(p).flip  = flip(p);
    
    
    if cbf(p)== 0
        cooked(:,p)=nan;
    else
        sim = gen_signals_vs_211024(parms(p),...
            timing_parms.t_delay,...
            timing_parms,...
            0,0);
        sim = abs(sim(1:Nframes))';
        
        cooked(:,p) = sim/norm(sim);
        raw(:,p) = raw(:,p) / norm(raw(:,p));
        
        if mod(p,1000)==0
            figure(777)
            
            plot(raw(:,p));
            hold on
            plot(cooked(:,p));
            hold off
            legend('raw','cooked')
            title(sprintf('vx %d of %d', p, Npix))
            drawnow
        end
        
    end
    
end

%%
error = mean(abs(raw - cooked), 1);
error = error ./ mean(abs(cooked), 1);
error(isnan(error)) = 0;
sub = raw-cooked;
sub(isnan(sub)) = 0;

save errors.mat raw cooked sub error

figure(778)

subplot(211)
imagesc(log(abs(sub(:,1:10:end))))
colorbar
title('Log Difference between model and data')
ylabel('time')
xlabel('voxels')

subplot(212)
lightbox(reshape(error,size(cbv)));
caxis([0 0.5])
title('residual')
colormap parula
