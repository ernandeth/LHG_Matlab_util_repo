function [dict, parms] = gen_t1dictionary_140328 (timing_parms)
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
%mtis0vals = linspace(700,1000,30);
mtis0vals = linspace(7000,12000,20);  % for JFN brain
mtis0vals = 1;

fvals = linspace(0.001, 0.02, 20);
fvals = 0.008;

cbvavals = linspace(0.01, 0.1, 10);
cbvavals = [0.02 0.5 0.9];

transvals = linspace(0.6, 2.5, 20);
transvals = 1;

kforvals = linspace(0.01, 0.5, 5);
kforvals = 0.3;
%kforvals = 0;

r1vals = 1./linspace(0.2, 4, 100);
%r1vals = 1/1.4;

betavals = (pi/180) * linspace(80, 100, 10);
betavals = (pi/180)*90;

Dispvals = linspace(0.01,10,10);
Dispvals = 10;

% % % % % REFERENCE % % % % %
% f =         0.01;
% mtis0 =     10e5 ;
% cbva =      0.04 ;
% transit=    1.0 ;
% kfor =      1e-3 ;
% r1tis =     1/1.4  ;

dofigs = 0;
doSub = 0;

for mtis0 = mtis0vals
    for transit = transvals
        for f = fvals
            for cbva = cbvavals
                for kfor = kforvals
                    for r1tis=r1vals;
                        for beta = betavals;
                            for Disp = Dispvals
                                
                            parms(n).mtis0 = mtis0;
                            parms(n).f = f;
                            parms(n).cbva = cbva;
                            parms(n).transit = transit;
                            parms(n).kfor = kfor;
                            parms(n).r1tis = r1tis;
                            parms(n).beta = beta;
                            parms(n).Disp = Disp;
                            
                            n = n+1;
                            %fprintf('\r .... %d',n);
                            end
                        end
                    end
                end
            end
        end
    end
end


entry = gen_signals_140328(parms(1) , timing_parms, 0,doSub);
Nobs = length(entry);
Ncombinations = n-1;
dict = zeros(Ncombinations, Nobs);

parms = parms(1:Ncombinations);

for n=1:Ncombinations;
%     fprintf('\rGenerating entry .... %d    of  %d   ',n, Ncombinations);
    entry = gen_signals_140328(parms(n) , timing_parms, dofigs,doSub);
    dict(n,:) = entry;
    parms(n);
end

if dofigs
    figure(2);
    imagesc(real(dict) );  title('dictionary entries (rows) Re')
    drawnow
end

save t1dictionary dict parms ;

end
