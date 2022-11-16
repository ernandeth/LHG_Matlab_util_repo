function [f2 transit2 cbv2] =  caslf_sim(SNR, useFileDictionary);

% Start out simple.... just generate a pulse seqeunce and time course
fprintf('\nTrue values for Single voxel case  ...\n')
f =         0.01;
mtis0 =     10e5 ;
cbva =      0.05 ;
transit=    1.2 ;   % this map will be different according to territory maps
kfor =      1e-2 ;      % MT parms still need to be verified
r1tis =     1/1.4  ;
L = 1;

if nargin==0
    SNR = 100;
    useFileDictionary = 0;
end

% variables in the dictionary entry with some default values
parms = struct( ...
    'mtis0', mtis0,...
    'f', f , ...  % perfusion in ml/s/g
    'cbva' ,  cbva , ...
    'transit', transit,...
    'kfor', kfor, ...
    'r1tis', r1tis);

% run the generator once to make a random acquisition sequence  and
% save the acquisition parameters for generating the dictionary and the
% observed signals
obs = gen_signals(parms, 0 , 1);
Nobs =length(obs);

%{
% Now let's deal with whole images ...

[gm h] = read_img('~hernan/matlab/SPM5/tpm/grey.nii');
[wm h] = read_img('~hernan/matlab/SPM5/tpm/white.nii');
[csf h] = read_img('~hernan/matlab/SPM5/tpm/csf.nii');
%
% % reduce their size to make is manageable:
% gm = gm(1:5:end, 1:5:end, 1:5:end);
% wm = wm(1:5:end, 1:5:end, 1:5:end);
% csf = csf(1:5:end, 1:5:end, 1:5:end);

% reduce their size to make is manageable:
gm = gm(1:5:end, 1:5:end, 36);
wm = wm(1:5:end, 1:5:end, 36);
csf = csf(1:5:end, 1:5:end, 36);

gm = gm/255;
wm = wm/255;
csf = csf/255;
L = length(csf(:));



%
% generate somw 'truth' maps
fprintf('\nGenerating synthetic parameter maps from GM, WM, CSF masks ...\n')
f =         gm * 0.008 + wm * 0.005 ;
mtis0 =     gm * 10e5  + wm * 10e5 +   csf*12e5 ;
cbva =      gm * 0.04 + wm * 0.02 ;
transit=    gm * 1.0 +  wm * 1.8;    % this map will be different according to territory maps
kfor =      gm * 1e-3 + wm * 1e-3  +   csf *1e-7;    % MT parms still need to be verified
r1tis =     gm /1.4  +  wm / 1.8  +   csf/1.4;

subplot(3,2,1) , lightbox(f);
subplot(3,2,2) , lightbox(mtis0);
subplot(3,2,3) , lightbox(cbva);
subplot(3,2,4) , lightbox(transit);
subplot(3,2,5) , lightbox(kfor);
subplot(3,2,6) , lightbox(r1tis);
%

drawnow
%}

useFileSignals=0;

if ~useFileSignals
    
    % allocate space
    parms = repmat(parms,[L 1]);
    % ... and  stick the physio parameters into the parms. structure
    for n=1:L
        fprintf('\r .... %d',n);
        parms(n).mtis0 = mtis0(n);
        parms(n).f = f(n);
        parms(n).cbva = cbva(n);
        parms(n).transit = transit(n);
        parms(n).kfor = kfor(n);
        parms(n).r1tis = r1tis(n);
    end
    fprintf('\n... parm structures filled for sythetic parameter maps');
    fprintf('\n genererating observed signals for those parameters.... \n');
    
    obssignals = zeros(L,Nobs);
    for n = 1:L
        %fprintf('\r pixel .... %d',n);
        obssignals(n,:) = gen_signals(parms(n), 1 , 0);
    end
    
    % add noise here:
    
    save observed obssignals
else
    load observed.mat
end

figure(5)
subplot(211)
plot(abs(obssignals)),title('signal magnitude')
subplot(212)
plot(angle(obssignals)); title('signal phase')

    
% Adding nois to the observation
% SNR = 400;
NoiseLevel = abs(mean(obssignals))/SNR
%NoiseLevel = 0
obssignals = obssignals + (NoiseLevel)*randn(size(obssignals));

figure(5)
subplot(211)
hold on, plot(abs(obssignals),'r'), hold off
subplot(212)
hold on, plot(angle(obssignals),'r'), hold off
drawnow

fprintf('\n... synthetic signals done');

% Dictionary Generation
% useFileDictionary=0;

if ~useFileDictionary
    fprintf('\nGenerating Dictionary ... \n')
    [dict, parms] = gen_dictionary;
    fprintf('\n ... done ')
else
    load dictionary.mat
end

fprintf('\n COmparing the observed signal at each voxel with each dictionary entry ... \n');
score = zeros(size(dict,1), 1);
P = size(obssignals,1);
for m=1:P
    
    fprintf('\r percent pixels .... %03f', 100*m/P);
    parfor n=1:size(dict,1);
        %fprintf('\r pixel .... %d   entry ... %d  ',m , n);
        cc = corrcoef(obssignals(m,:), dict(n,:));
        score(n) = abs(cc(1,2));
        %score(n) = obssignals(m,:) * dict(n,:)';
        %
    end
    best = max(score(:));
    
    % fprintf('\ngrabbing the best fits...');
    if ~isnan(best)
        bestind = find(score==best);
        if length(bestind)>1
            fprintf('\nMore than one entry fits the maximum value ... put in Nan')
            bestind=nan;
            bestparms(m) = nan;
        else
            bestparms(m) = bestind;
        end
    else
        bestparms(m) = nan;
    end
end


save bestparmsinds bestparms score
fprintf('\n .... done');

for n=1:length(bestparms)
    if ~isnan(bestparms(n))
        f2(n) = parms(bestparms(n)).f;
        transit2(n) = parms(bestparms(n)).transit;
		cbv2(n) = parms(bestparms(n)).cbva;
    else
        f2(n) = 0;
        transit2(n) = 0;
		cbvs2(n) = 0;
    end
end

save truth_images f mtis0 transit cbva kfor r1tis

parms(bestparms)
figure(3); subplot(211)
plot(score), title('Correlation with dictionary entries')
subplot(212)
hist(score,100);
%f2 = reshape(f2,size(csf));
%lightbox(f2),[0 0.0078],4;
% pause

return


function [dict, parms] = gen_dictionary
%%
% this is the function that generates a dictionary
%
% First we need to create the dictionary of signals
% this requires creating a model for how each entry's signal
% should look.  Use Bloch simulation for this

% variables in the dictionary entry with some default values
parms = struct( ...
    'mtis0', 0,...
    'f', 0 , ...  % perfusion in ml/s/g
    'cbva' ,  0 , ...
    'transit', 0,...
    'kfor', 0, ...
    'r1tis', 0);

% the parameters will be an array of 5 cells

%fprintf('\nAllocating space for dictionary ... \n');
L = 1e5;
parms = repmat(parms,[L 1]);

n = 1;

% mtis0vals = 1e5*[8 10 20];
% fvals = [0.003:0.001:0.01];
% cbvavals = [0.01:0.01:0.05];
% transvals = [0.8:0.1:1.5];
% kforvals = [0 0.01 0.05];
% r1vals = 1./[1: 0.1: 2];

% this is for the gray matter only
% mtis0vals = 1e5*linspace(6,8,5);
mtis0vals = 10e5;
fvals = linspace(0.001, 0.015, 50);
% fvals = 0.01;
cbvavals = linspace(0.01, 0.1, 10);
%cbvavals = 0.05;
transvals = linspace(0.6, 2, 50);
% transvals = 1.0;
% kforvals = linspace(0,0.01,5);
kforvals = 1e-2;
% r1vals = 1./linspace(1,1.4,5);
r1vals = 1/1.4;

% for reference:  this is how I made the 'truth' maps above.
% fprintf('\nTrue values for Single voxel case  ...\n')
% f =         0.01;
% mtis0 =     10e5 ;
% cbva =      0.04 ;
% transit=    1.0 ;   % this map will be different according to territory maps
% kfor =      1e-3 ;      % MT parms still need to be verified
% r1tis =     1/1.4  ;


for mtis0 = mtis0vals   
    for transit = transvals
        for f = fvals
            for cbva = cbvavals
                for kfor = kforvals
                    for r1tis=r1vals;
                        parms(n).mtis0 = mtis0;
                        parms(n).f = f;
                        parms(n).cbva = cbva;
                        parms(n).transit = transit;
                        parms(n).kfor = kfor;
                        parms(n).r1tis = r1tis;
                        n = n+1;
                        fprintf('\r .... %d',n);
                    end
                end
            end
        end
    end
end

entry = gen_signals(parms(1) , 1, 0);
Nobs = length(entry);
Ncombinations = n-1;
dict = zeros(Ncombinations, Nobs);
dofigs = 0;
parms = parms(1:Ncombinations);

for n=1:Ncombinations;
    fprintf('\rGenerating entry .... %d    of  %d   ',n, Ncombinations);
    entry = gen_signals(parms(n) , 1, dofigs);
    dict(n,:) = entry;
end

figure(2);
subplot(211);
imagesc(real(dict) );  title('dictionary entries (rows) Re')
subplot(212);
%imagesc(abs(corr(dict'))); title('correlations among entries')
drawnow
save dictionary dict parms ;


return



function obs = gen_signals(parms, getAQfromfile, dofigs)
%%
% set up parms for pulse sequence and known constants
% ... time units are in seconds

duration = 500;
dt = 0.01;
Disp = 0.1/dt;
r1a = 1/1.6;
TR = 4;

% get parms from the input structure
mtis0 = parms.mtis0;
f = parms.f;
cbva = parms.cbva;
transit = parms.transit;
kfor = parms.kfor;
r1tis = parms.r1tis;

% seq timings, the sequence cycle (TR) is 4 seconds
startL = [2:TR*2:duration-15];
startC = startL + TR;
aqtime = [startL startC];
aqtime = sort(aqtime)-0.4;

N = length(startL);

doBS = 1;

if ~getAQfromfile
    fprintf('\nGenerating a new random acquisition sequence ...');
    option=3
    
    mt_aqs = [];
    
    switch(option)
        case 1
            % random label times and durations, but fixed BS pulse timing
            % this one does worst
            label_dur = 3*(rand(size(startL)));
            control_dur = 3*(rand(size(startC)));
            for n=1:length(label_dur)
                label_start(n) = (3 - label_dur(n))*rand(1);
                control_start(n) = (3 - control_dur(n))*rand(1);
            end
            
            bstime0 = zeros(size(aqtime)) + 0.02;
            bstime1 = 1.9 * ones(size(aqtime));
            bstime2 = 3.0 * ones(size(aqtime));
            
            
        case 2
            % or a sequential increase in labeling times
            % this one does best
            label_dur = linspace(0,3,length(startL));
            control_dur = linspace(0,3,length(startC));
            label_start = 3 - label_dur;
            control_start = 3 - control_dur;
            
            bstime0 = zeros(size(aqtime));
            bstime1 = 2.1 * ones(size(aqtime));
            bstime2 = 3.0 * ones(size(aqtime));
            
        case 3
            % random label times and durations, but fixed BS pulse timing
            % label and control are matched
            label_dur = 3*(rand(size(startL)));
            
            for n=1:length(label_dur)
                label_start(n) = (3 - label_dur(n))*rand(1);
            end
            
            control_dur = label_dur;
            control_start = label_start;
            
            bstime0 = zeros(size(aqtime)) + 0.02;
            bstime1 = 1.9 * ones(size(aqtime));
            bstime2 = 3.0 * ones(size(aqtime));
             
       
    end
    
    % in these Aq periods, the label will saturate spins, not
    % invert.
    mt_aqs = randperm(N);
    mt_aqs = mt_aqs(1:end/4);
    
    
    
    save AQparms label_start control_start label_dur control_dur bstime1 bstime2 bstime0 aqtime mt_aqs
else
    %fprintf('\rReading acquisition sequence from file...');
    load AQparms
end

labelfun = zeros(duration/dt,1);
bsfun = zeros(duration/dt,1);
aqfun = zeros(duration/dt,1);
t = linspace(0,duration,length(labelfun));

% labeling function can have 3 values: off(0), label (1), control(-1)

for n=1:length(startL)
    
    ind1 = 2 + (2*n-2)*TR + label_start(n);
    ind2 = 2 + (2*n-2)*TR + label_start(n) + label_dur(n);
    ind3 = 2 + (2*n-2)*TR + bstime0(n);
    ind4 = 2 + (2*n-2)*TR + bstime1(n);
    ind5 = 2 + (2*n-2)*TR + bstime2(n);
    
    labelfun( round(ind1/dt): round(ind2/dt) ) = 1;
    bsfun(round([ind3 ind4 ind5]/dt)) = 1;
    
    % make the label be a saturation during some cases.
    if ~isempty(find(mt_aqs==n))
        labelfun( round(ind1/dt): round(ind2/dt) ) = 0.5;
    end
    
    ind1 = 2 + (2*n-1)*TR + control_start(n);
    ind2 = 2 + (2*n-1)*TR + control_start(n) + control_dur(n);
    ind3 = 2 + (2*n-1)*TR + bstime0(n);
    ind4 = 2 + (2*n-1)*TR + bstime1(n);
    ind5 = 2 + (2*n-1)*TR + bstime2(n);
    
    labelfun( round(ind1/dt): round(ind2/dt) ) = -1;
    bsfun(round([ind3 ind4 ind5]/dt)) = 1;
    
    % make the label be a saturation during some cases.
    if ~isempty(find(mt_aqs==n))
        labelfun( round(ind1/dt): round(ind2/dt) ) = 0.5;
    end
    
end


bsfun(2) = 1;
aqfun(round(aqtime/dt)) = 1;
aqs = find(aqfun);

%{
area(t,labelfun)
hold on
plot(t,bsfun,'r')
hold off
%}

%%% remove background suppression  ****
if ~doBS
    bsfun(:) = 0;
end
%%%%%%%%%%%%%%%%

% form an arterial dispersion kernel (only need to do this once)
ttmp = linspace(0,5, 5/dt)';
art_kernel = (ttmp-transit).* exp(-Disp *(ttmp-transit));
art_kernel(art_kernel<0)=0;
art_kernel = art_kernel / sum(art_kernel);

% art_kernel = art_kernel .* exp(- r1a*ttmp);
% plot(ttmp, art_kernel)

% calculating the arterial signals
% The input function into the voxel is the label after t1 decau and
% dispersion
inp = ones(size(labelfun));
inp(labelfun==1) = 1 - 2*0.8* exp(-transit*r1a);
inp(labelfun==0.5) = 1 - exp(-transit*r1a);  % in this case we get arterial spin saturation
tmp_inp = inp;
%plot(inp)

% pad the input function with ones
inp = [ones(size(ttmp)); inp];

% convolve with dipersion kernel
% in order to obtain the arterial signal
% assume that the arterial volume is 100% for input function
ma = mtis0 * conv(inp, art_kernel);

% remove padding from both ends
ma = ma(length(ttmp):end);
ma = ma(1:length(labelfun));
ma = ma * mtis0/max(ma);

%{
plot(t,ma)
hold on
plot(t, ma, 'r')
hold off
%}

% now the fun part: calculate the tissue function
% make a matgnetization transfer function

mtr = kfor * ~~(labelfun);  % we always get the same RF power and MT whether there is label or not.
m = mtis0 * ones(size(t));

for n=2:length(t)
    % the modified Blocj equation has
    % t1 decay,  magnetization transfer, inflow, outflow
    dm = (mtis0 - m(n-1))*r1tis - mtr(n-1)*m(n-1) + f*ma(n-1) - f*m(n-1)/0.9;
    m(n) = m(n-1) + dm*dt;
    
    if bsfun(n)==1
        m(n) = -m(n);
        ma(n) = -ma(n);
    end
    
end

% adjust ma to reflect the blood volume, not just for a 100% voxel
ma = cbva*ma';
%m = m + ma;

obs = m(aqs) + ma(aqs);
% modification:  signals are complex.   Assume flip angle beta is pi/2
% Arterial and tissue compartments will acquire different amounts of
% phase during readout
beta = pi/2;
m_phi = 0.1;
ma_phi = 2;
obs = m(aqs)*sin(beta)*exp(-i*m_phi) + ma(aqs)*sin(beta)*exp(-i*ma_phi);
obs = obs(3:end);
%obs = ma(aqs);

if dofigs
    figure(1)
    subplot(411)
    area(t,labelfun)
    hold on
    plot(t,tmp_inp)
    plot(t,bsfun,'r')
    plot(t(aqs),bsfun(aqs),'g*')
    title('label/control/BS pulses')
    hold off
    
    subplot(412)

    hold on
    plot(t,ma,'r')
    plot(t(aqs),ma(aqs),'g*')
    hold off
    title('arterial magnetization')
    
    subplot(413)
    plot(t,m)
    hold on
    plot(t(aqs),m(aqs),'g*')
    hold off
    title('tissue magnetization')
    
    subplot(414)
    plot(angle(obs))

    title('phase of aquired signal')
    
    
    figure(5)
    subplot(211)
    plot(abs(obs)); title('signal magnitude')
    subplot(212)
    plot(angle(obs)); title('signal phase')
    
    
    drawnow
end



return
