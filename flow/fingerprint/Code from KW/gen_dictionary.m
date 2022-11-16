function [dict, parms] = gen_dictionary
% Use Bloch simulation for this

parms = struct( ...
    'mtis0', 0,...
    'f', 0 , ...  
    'cbva' ,  0 , ...
    'transit', 0,...
    'kfor', 0, ...
    'r1tis', 0);

L = 1e5;
parms = repmat(parms,[L 1]);

n = 1;
% biiiig dictionary: 
% mtis0vals = 1;
% fvals = 0.001:0.001:0.015;
% cbvavals = 0:0.1:1;
% transvals = 0.5:0.1:2.0;
% kforvals = 1e-4;
% r1vals = 1./(0.8:0.05:1.6);

% Huge dictionary created at 730pm
% mtis0vals = 1;%0e5;
% % fvals=0.001;
% % fvals = 0.001:0.001:0.02;
% fvals = 0.001:0.001:0.02;
% cbvavals = 0:0.1:1;
% % cbvavals = linspace(0.01, 0.1, 10);
% % transvals = 0.5:0.05:2.5;%linspace(0.6, 2, 10);
% transvals = 0.3:0.05:2.0;%linspace(0.6, 2, 10);
% % transvals = 1.5;
% % kforvals = linspace(0,0.01,5);
% kforvals = 1e-4;
% % r1vals = 1./(0.4:0.05:1.6);%1./linspace(.2,2,50);
% r1vals = 1./(0.4:0.05:1.8);%1./linspace(.2,2,50);
% r1vals = 1/1.4;

mtis0vals = 1;%0e5;
% fvals=0.001;
% fvals = 0.001:0.001:0.02;
fvals = 0.001:0.0025:0.02;
cbvavals = 0:0.2:1;
% cbvavals = linspace(0.01, 0.1, 10);
% transvals = 0.5:0.05:2.5;%linspace(0.6, 2, 10);
transvals = 0.3:0.1:2.0;%linspace(0.6, 2, 10);
% transvals = 1.5;
% kforvals = linspace(0,0.01,5);
kforvals = 1e-4;
% r1vals = 1./(0.4:0.05:1.6);%1./linspace(.2,2,50);
r1vals = 1./(0.4:0.1:1.8);%1./linspace(.2,2,50);

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

entry = gen_signals_110613(parms(1) , 1, 0);
Nobs = length(entry);
Ncombinations = n-1;
dict = zeros(Ncombinations, Nobs);
dofigs = 0;
parms = parms(1:Ncombinations);

matlabpool(4)
parfor n=1:Ncombinations;
    fprintf('\rGenerating entry .... %d    of  %d   ',n, Ncombinations);
    entry = gen_signals_110613(parms(n) , 1, dofigs);
    dict(n,:) = entry;
end
matlabpool close;
figure(2);
% subplot(211);
imagesc(real(dict) );  title('dictionary entries (rows) Re')
% subplot(212);
% drawnow
save dictionary dict parms ;

end
