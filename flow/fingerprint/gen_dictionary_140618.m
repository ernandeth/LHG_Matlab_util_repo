function [dict, parms] = gen_dictionary_140618 (timing_parms)
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

mtis0vals = 1;

fvals = linspace(0, 100,4 )/6000;
%fvals = 40/6000;
%fvals = 0;

cbvavals = [linspace(0.001, 0.01,3) ];
%cbvavals = 0.01;
%cbvavals = 0;

transvals = linspace(0.5, 3, 3);
%transvals = 1.2;

%kforvals = linspace(0.0,1,5);
kforvals = 0.02;

r1vals = 1./linspace(0.5,3,10);
%r1vals = 1./[0.7:0.1:2  3];
r1vals = 1./[0.9 1.4 3];
r1vals = 1/1.4;

betavals =linspace(20,60,10) * pi/180;
betavals = deg2rad(70 );

%Dispvals = linspace(0.01,40,10);
%Dispvals = [10 20 30 40 ];
Dispvals = [50];

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

entry = gen_signals_140618(parms(1) , timing_parms, dofigs, doSub);
Nobs = length(entry);
Ncombinations = n-1;
dict = zeros(Ncombinations, Nobs);
parms = parms(1:Ncombinations);

%matlabpool open 4
parfor n=1:Ncombinations;
%    fprintf('\rGenerating entry .... %d    of  %d   ',n, Ncombinations);
    entry = gen_signals_140618(parms(n) , timing_parms, dofigs, doSub);
    dict(n,:) = entry;
    parms(n);
end
%matlabpool close

if dofigs
    figure(2);
    imagesc(real(dict) );  title('dictionary entries (rows)')
    drawnow
end

save dictionary dict parms ;

end
