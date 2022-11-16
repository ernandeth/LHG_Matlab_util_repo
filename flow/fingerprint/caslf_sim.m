function [f2 transit2 cbv2] =  caslf_sim(SNR, useFileDictionary);

% fprintf('\nTrue values for Single voxel case  ...\n')
f =         0.0096;
mtis0 =     10e2 ;
cbva =      0.05 ;
transit=    1.2 ;   % this map will be different according to territory maps
kfor =      1e-2 ;      % MT parms still need to be verified
r1tis =     1/1.4  ;
L = 1;
doFigs=0;


% default values
parms = struct( ...
    'mtis0', mtis0,...
    'f', f , ...
    'cbva' ,  cbva , ...
    'transit', transit,...
    'kfor', kfor, ...
    'r1tis', r1tis);

if nargin==0
    SNR=500;
    useFileDictionary=1;
end

useFileSignals=0;
obs = gen_signals(parms, 0 , 1);
Nobs =length(obs);

if ~useFileSignals
    
    parms = repmat(parms,[L 1]);
    
    for n=1:L
        %         fprintf('\r .... %d',n);
        parms(n).mtis0 = mtis0(n);
        parms(n).f = f(n);
        parms(n).cbva = cbva(n);
        parms(n).transit = transit(n);
        parms(n).kfor = kfor(n);
        parms(n).r1tis = r1tis(n);
    end
    %     fprintf('\n... parm structures filled for sythetic parameter maps');
    %     fprintf('\n genererating observed signals for those parameters.... \n');
    
    obssignals = zeros(L,Nobs);
    for n = 1:L
        obssignals(n,:) = gen_signals(parms(n), 1 , 1);
    end
    
    
    save observed obssignals
else
    load observed.mat
end

% Adding noise to the observation
NoiseLevel = abs(mean(obssignals))./SNR;
obssignals = obssignals + (NoiseLevel).*randn(size(obssignals));

if doFigs
    figure(5)
    subplot(211)
    hold on, plot(abs(obssignals),'r'), hold off
    subplot(212)
    hold on, plot(angle(obssignals),'r'), hold off
    drawnow
end

% fprintf('\n... synthetic signals done');

bestparms = search_dictionary(obssignals, useFileDictionary, doFigs)

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

end






