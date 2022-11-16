function [dict, parms] = gen_dictionary
% Use Bloch simulation for this
doSub=0;
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

% this is for the gray matter only
mtis0vals = 10e2;

fvals = linspace(0.001, 0.02, 50);

%cbvavals = 0.05;
cbvavals = linspace(0.001, 0.1, 10);

%transvals =1.2;
transvals = linspace(0.8, 1.8, 50);

% kforvals = linspace(0,0.01,5);
kforvals = 1e-2;

% r1vals = 1./linspace(1,1.4,5);
r1vals = 1/1.4;

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
                        %fprintf('\r .... %d',n);
                    end
                end
            end
        end
    end
end


entry = gen_signals(parms(1) , doSub, 0);
Nobs = length(entry);
Ncombinations = n-1;
dict = zeros(Ncombinations, Nobs);
dofigs = 0;
parms = parms(1:Ncombinations);

for n=1:Ncombinations;
    fprintf('\rGenerating entry .... %d    of  %d   ',n, Ncombinations);
    entry = gen_signals(parms(n) , doSub, dofigs);
    dict(n,:) = entry;
    parms(n)
end

figure(2);
subplot(211);
imagesc(real(dict) );  title('dictionary entries (rows) Re')
subplot(212);
drawnow
save dictionary dict parms ;

end
