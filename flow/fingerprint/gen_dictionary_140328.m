function [dict, parms] = gen_dictionary_140328 (timing_parms)
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
%mtis0vals = linspace(7000,12000,20);  % for JFN brain
mtis0vals = 1;

fvals = linspace(0, 300, 100)/6000;
% fvals = [0.005 0.01];

cbvavals = linspace(0.01, 0.5, 5);
%cbvavals = 0.02;

transvals = linspace(0.2, 1.2, 10);
%transvals = 0.8;

kforvals = linspace(0.15,0.35,20);
kforvals = 0.2;

r1vals = 1./linspace(0.5,2,5);
r1vals = 1./[0.9 1.8 3];
r1vals = 1/1.7;

%betavals = (pi/180) * linspace(15,90,20);
betavals = (pi/180)*[50 60 70 80 ];
betavals = 90 * pi/180;

Dispvals = linspace(0.01,40,10);
Dispvals = [10 20 30 40 ];
Dispvals = 13;

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

doSub = 0;
dofigs = 0;

entry = gen_signals_140328(parms(1) , timing_parms, dofigs, doSub);
Nobs = length(entry);
Ncombinations = n-1;
dict = zeros(Ncombinations, Nobs);
parms = parms(1:Ncombinations);

for n=1:Ncombinations;
%     fprintf('\rGenerating entry .... %d    of  %d   ',n, Ncombinations);
    entry = gen_signals_140328(parms(n) , timing_parms, dofigs, doSub);
    dict(n,:) = entry;
    parms(n);
end

if dofigs
    figure(2);
    imagesc(real(dict) );  title('dictionary entries (rows)')
    drawnow
end

save dictionary dict parms ;

end
