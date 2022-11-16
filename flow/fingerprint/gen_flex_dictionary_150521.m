function [dict, parms] = gen_flex_dictionary_150521 (timing_parms, fixed_parms)
% function [dict, parms] = gen_flex_dictionary_150521 (timing_parms, [fixed_parms])

% % % % % REFERENCE % % % % %
% f =         0.01;
% mtis0 =     10e5 ;
% cbva =      0.04 ;
% bat=    1.0 ;
% kfor =      1e-3 ;
% r1tis =     1/1.4  ;

% allocate space:
parms = struct( ...
    'mtis0', 0,...
    'f', 0 , ...
    'cbva' ,  0 , ...
    'bat', 0,...
    'kfor', 0, ...
    'r1tis', 0,...
    'flip',0, ...
    'Disp',0);


n = 1;

% default search space to construct dictionary
mtis0vals = 1;
fvals = linspace(0, 100, 20) / 6000;
cbvavals = [linspace(0.001, 0.01, 20) 1];
batvals = linspace(0.5, 3, 20);
kforvals = linspace(0,0.1,20);
r1vals = linspace(0.3,2,20);
flipvals = linspace(20,60,20) * pi/180;
Dispvals = linspace(1,50,20);

if nargin==2
    
    if isfield(fixed_parms, 'mtis0'),  mtis0vals = fixed_parms.mtis0;  end
    
    if isfield(fixed_parms, 'f') ,  fvals = fixed_parms.f;  end
    
    if isfield(fixed_parms, 'cbva') ,  cbvavals = fixed_parms.cbva;  end
    
    if isfield(fixed_parms, 'bat') ,  batvals = fixed_parms.bat;  end

    if isfield(fixed_parms, 'bat2') ,  bat2vals = fixed_parms.bat2;  end

    if isfield(fixed_parms, 'kfor'),  kforvals = fixed_parms.kfor;  end
    
    if isfield(fixed_parms, 'r1tis'),  r1vals = fixed_parms.r1tis;  end
    
    if isfield(fixed_parms, 'flip'),  flipvals = fixed_parms.flip;  end
    
    if isfield(fixed_parms, 'Disp'),  Dispvals = fixed_parms.Disp;  end
    
end

Ncombinations = length(Dispvals) * length(flipvals) * length(r1vals) * ...
    length(kforvals) * length(batvals)* length(bat2vals) * length(cbvavals) * length(fvals) * length(mtis0vals) ;

%keyboard

parms = repmat(parms,[Ncombinations 1]);

for mtis0 = mtis0vals
    for bat = batvals
        for bat2 = bat2vals
            
            for f = fvals
                for cbva = cbvavals
                    for kfor = kforvals
                        for r1tis=r1vals;
                            for flip = flipvals;
                                for Disp = Dispvals
                                    
                                    parms(n).mtis0 = mtis0;
                                    parms(n).f = f;
                                    parms(n).cbva = cbva;
                                    parms(n).bat = bat;
                                    parms(n).bat2 = bat2;
                                    parms(n).kfor = kfor;
                                    parms(n).r1tis = r1tis;
                                    parms(n).flip = flip;
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
end

doSub = 0;
dofigs = 0;

% entry = gen_signals_150507(parms(1) , timing_parms, dofigs, doSub);
%entry = gen_signals_150521(parms(1) , timing_parms, dofigs, doSub);
entry = gen_signals_160426(parms(1) , timing_parms, dofigs, doSub);

Nobs = length(entry);
Ncombinations = n-1;
dict = zeros(Ncombinations, Nobs, 'single');
parms = parms(1:Ncombinations);

%matlabpool open 4
parfor n=1:Ncombinations;
    %    fprintf('\rGenerating entry .... %d    of  %d   ',n, Ncombinations);
    %entry = gen_signals_150507(parms(n) , timing_parms, dofigs, doSub);
    % entry = gen_signals_150521(parms(n) , timing_parms, dofigs, doSub);
    entry = gen_signals_160426(parms(n) , timing_parms, dofigs, doSub);

    
    % Normalize all the entries of the dictionary
    tmp = entry - mean(entry);
    tmp = tmp / norm(tmp);
    entry = tmp;
    %
    dict(n,:) = entry;
    parms(n);
end
%matlabpool close


if dofigs
    figure(2);
    imagesc(real(dict) );  title('dictionary entries (rows)')
    drawnow
end

%save dictionary dict parms ;

end
