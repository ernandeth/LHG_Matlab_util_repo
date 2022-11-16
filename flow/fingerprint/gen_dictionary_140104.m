function [dict, parms] = gen_dictionary_1401014
% make a dictionary for fitting T1, flipangle and M0

parms = struct( ...
    'mtis0', 0,...
    'f', 0 , ...
    'cbva' ,  0 , ...
    'transit', 0,...
    'kfor', 0, ...
    'r1tis', 0,...
    'beta',0);

L = 1e5;
parms = repmat(parms,[L 1]);

n = 1;

% this is for the gray matter only
mtis0vals = linspace(700,1000,30);
mtis0vals = linspace(7000,12000,20);  % for JFN brain

%fvals = linspace(0.001, 0.02, 50);
fvals = 0;

%cbvavals = linspace(0.001, 0.1, 10);
cbvavals = 0;

% transvals = linspace(0.8, 1.8, 50);
transvals = 0;

kforvals = linspace(0,0.20,10);
% kforvals = 1e-2;
% kforvals = 0;

r1vals = 1./linspace(0.9,1.4,10);
%r1vals = 1/1.4;

betavals = (pi/180) * linspace(15,90,20);
% betavals = (pi/180)*90;

% % % % % REFERENCE % % % % %
% f =         0.01;
% mtis0 =     10e5 ;
% cbva =      0.04 ;
% transit=    1.0 ;
% kfor =      1e-3 ;
% r1tis =     1/1.4  ;


for mtis0 = mtis0vals
    for transit = transvals
        for f = fvals
            for cbva = cbvavals
                for kfor = kforvals
                    for r1tis=r1vals;
                        for beta = betavals;
                            parms(n).mtis0 = mtis0;
                            parms(n).f = f;
                            parms(n).cbva = cbva;
                            parms(n).transit = transit;
                            parms(n).kfor = kfor;
                            parms(n).r1tis = r1tis;
                            parms(n).beta = beta;
                            n = n+1;
                            %fprintf('\r .... %d',n);
                        end
                    end
                end
            end
        end
    end
end

entry = gen_signals_140104(parms(1) , 1, 0);
Nobs = length(entry);
Ncombinations = n-1;
dict = zeros(Ncombinations, Nobs);
dofigs = 0;
parms = parms(1:Ncombinations);

for n=1:Ncombinations;
    fprintf('\rGenerating entry .... %d    of  %d   ',n, Ncombinations);
    entry = gen_signals_140104(parms(n) , 1, dofigs);
    dict(n,:) = entry;
    parms(n);
end

figure(2);
subplot(211);
imagesc(real(dict) );  title('dictionary entries (rows) Re')
subplot(212);
drawnow
save dictionary dict parms ;

end
