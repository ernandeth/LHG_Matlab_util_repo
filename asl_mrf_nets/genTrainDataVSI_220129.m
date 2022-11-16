% Generates ASL Training Data for a machine learning based estimation
% technique. The parameters are uniformly random over a range, 
%
% evolved from AL code for VSI training
% LHG: jan 2022
% modiifcations for new VSASL model.  Key difference:  reduced number of
% parameters.  (No kfor, bat2 ... alphas)
%
% All the labeling efficiencies are calculated from T2 of the 
% blood and tissue, and the type of pulse. 
% - Note that each  type VS pulse has an effective TE

useFileSchedule = 1;

% number of training + testing time-series and the length of time series;
nData = 6.0*10^6;
nData = 2.0*10^6; % New populations: 2/17/22

% nData = 50;  % small data set to make sure things work right


if useFileSchedule
    % read in the timing paramaters from text files
    formatSpec = '%f';
    fileID = fopen('AL_schedule_model_220130/t_delays.txt','r');
    timing_parms.t_delay = fscanf(fileID,formatSpec);
    fileID = fopen('AL_schedule_model_220130/t_adjusts.txt','r');
    timing_parms.t_adjusts = fscanf(fileID,formatSpec);
    fileID = fopen('AL_schedule_model_220130/doArtSuppression.txt','r');
    timing_parms.doArtSup =  fscanf(fileID,formatSpec);
    fileID = fopen('AL_schedule_model_220130/isVelocitySelective.txt','r');
    timing_parms.labelcontrol =  fscanf(fileID,formatSpec);
    
    % some preset constant values:
    timing_parms.label_type     = 'FTVSI-sinc'; %'FTVSI-sinc'; % 'BIR8inv'; % 'BIR8'
    timing_parms.readout_type   = 'FSE'; %'GRE'  % 3/16/22
    Nframes                     = length(timing_parms.t_delay);
    timing_parms.t_tag          = 0.0*ones(Nframes,1) ;
    timing_parms.ArtSup_delay   = 0.100*ones(Nframes,1) ;
    %timing_parms.t_aq           = 0.0329*18*ones(Nframes,1) ;
    timing_parms.t_aq           = (0.031*18 +0.016)*ones(Nframes,1) ; %3/16/22
    timing_parms.order          = 1;

else
    % make up our own schedule on the fly
    Nframes = 150;
    rp  = randperm(Nframes); rp=rp(1:end/2);
    rp2 = randperm(Nframes); rp2=rp2(1:end/2);
    
    timing_parms.labelcontrol           = ones(Nframes,1);
    timing_parms.labelcontrol(rp)  = 0;
    timing_parms.doArtSup               = ones(Nframes,1);
    timing_parms.doArtSup(rp2)      = 0;
    
    timing_parms.t_adjusts              = 0.1 *ones(Nframes,1);
    timing_parms.t_delay                = 0.2 + 1.3*abs(chirp(linspace(0.75,0,Nframes), 0, 1,3.5))';
    timing_parms.t_delay(2:2:end)       = timing_parms.t_delay(1:2:end);
    %timing_parms.t_delay(:)       = 1;
    
    timing_parms.t_tag                  = 0.10*ones(Nframes,1) ;
    timing_parms.ArtSup_delay           = 0.100*ones(Nframes,1) ; % verified on 3/2/22
    timing_parms.t_aq                   = 0.0329*18*ones(Nframes,1) ;
    timing_parms.order                  = 1;
    timing_parms.label_type             = 'FTVSI-sinc'; %'FTVSI-sinc'; % 'BIR8inv'; % 'BIR8'
    timing_parms.readout_type           = 'FSE';

end

%2/22/22 : let's see if I'm wrong about the order of labels
% timing_parms.order          = 2;
%%%%

% synthetic training data matrix
Y = single(zeros(nData,Nframes)); 

% train-test split ratio
split = 0.9;

%
delays = timing_parms.t_delay;

seqDuration = sum(timing_parms.t_delay)+ ...
    sum(timing_parms.t_tag) + ...
    sum(timing_parms.t_adjusts) + ...
    sum(timing_parms.ArtSup_delay)+...
    sum(timing_parms.t_aq)
%
%figures and subtraction
doFigs = 0;
doSub = 0;

% generate Data : single run
tic

parms.mtis0 =   1 ;
parms.Disp =    40;
parms.r1blood = 1/1.7;

parms.f =       0.005;
parms.cbva =    0.005;
parms.bat =     0.01;
parms.r1tis =   1/1.2;
parms.r2tis =   1/0.050;
parms.flip =    deg2rad(30);

test_signal = abs(single(  ...
    gen_signals_vs_211024(parms,...
    delays, ...
    timing_parms,...
    doFigs,doSub)));

plot(test_signal)
drawnow
%%
% Now generate distribution of tissue parameters
%{
order:
    parms.f =       p_array(1);
    parms.cbva =    p_array(2);
    parms.bat =     p_array(3);
    parms.r1tis =   p_array(4);
    parms.r2tis =   p_array(5);
    parms.flip =    p_array(6);

%}
f       = rand(1,nData) * 0.02;
cbv     = rand(1,nData) * 0.05;
bat_tis = 0.05 + rand(1,nData) * 1.2;
%2/17/22 : New population :
%r1tis   = 0.3 + 2*rand(1,nData) ;
%r2tis   =  10 + 40*rand(1,nData) ;
%flip    =  pi/20 + (pi/2)*rand(1,nData);
r1tis   =  0.25 + 3*rand(1,nData) ;
r2tis   =  5 + 45*rand(1,nData) ;
flip    =  0.1 + (pi/4)*rand(1,nData);

% bat_tis = 0.05 + rand(1,nData) * 0.5;  % 3/1/22

theta = [
    f;
    cbv;
    bat_tis;
    r1tis;
    r2tis;
    flip]; 

theta = single(theta)';


delays = timing_parms.t_delay;



%%
fprintf('\nBegin generating %d time courses', nData);
%parpool('local', 'SpmdEnabled', false)
%pp =  parpool(24)
doFigs = 0;

parfor idx = 1:nData
    
    parms = make_parm_struct(theta(idx,:));
    
    Y(idx,:) = abs(single(  ...
        gen_signals_vs_211024(...
        parms,...
        delays,...
        timing_parms,...
        doFigs,doSub)));
        
    if mod(idx,1e3)==0
        fprintf('\ngenerating time course number ... %d', idx);
    end
end

fprintf('\n ....Finished generating %d time courses', nData);

toc

%% split into training and testing data
tic
lim = round(split*nData);
trainData = Y(1:lim,:);
trainTgts = theta(1:lim,:);

testData = Y(lim+1:end,:);
testTgts = theta(lim+1:end,:);
%%

!mkdir AL_schedule_model_220130/
save -v7.3 AL_schedule_model_220130/trainData_4.mat trainData
save -v7.3 AL_schedule_model_220130/trainTgts_4.mat trainTgts
save -v7.3 AL_schedule_model_220130/testData_4.mat testData
save -v7.3 AL_schedule_model_220130/testTgts_4.mat testTgts
toc

%%
function parms = make_parm_struct(p_array)

    parms.f =           p_array(1);
    parms.cbva =        p_array(2);
    parms.bat =         p_array(3);
    parms.r1tis =       p_array(4);
    parms.r2tis =       p_array(5);
    parms.flip =        p_array(6);
    parms.mtis0 =       1 ;
    parms.Disp =        40;
    parms.r1blood =     1/1.7;

    
end
