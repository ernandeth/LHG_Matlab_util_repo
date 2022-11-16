function [dict, parms] = gen_flex_dictionary_160525 (timing_parms, fixed_parms)
% function [dict, parms] = gen_flex_dictionary_160525 (timing_parms, [fixed_parms])

% % % % % REFERENCE % % % % %
% f =         0.01;
% mtis0 =     10e5 ;
% cbva =      0.04 ;
% bat=    1.0 ;
% bat2=    1.0 ;
% kfor =      1e-3 ;
% r1tis =     1/1.4  ;




n = 1;

% default search space to construct dictionary
mtis0vals = 1;
fvals = linspace(0, 100, 20) / 6000;
cbvavals = [linspace(0.001, 0.01, 20) 1];
batvals = linspace(0.5, 3, 20);
bat2vals = linspace(0.5, 3, 20);
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

% allocate space:
parms = struct( ...
    'mtis0', single(zeros(Ncombinations,1)),...
    'f', single(zeros(Ncombinations,1)) , ...
    'cbva' ,  single(zeros(Ncombinations,1)) , ...
    'bat', single(zeros(Ncombinations,1)),...
    'bat2', single(zeros(Ncombinations,1)),...
    'kfor', single(zeros(Ncombinations,1)), ...
    'r1tis', single(zeros(Ncombinations,1)),...
    'flip',single(zeros(Ncombinations,1)), ...
    'Disp',single(zeros(Ncombinations,1)) );
%keyboard

%parms = repmat(parms,[Ncombinations 1]);

for mtis0 = mtis0vals
    for bat = batvals
        for bat2 = bat2vals          
            for f = fvals
                for cbva = cbvavals
                    for kfor = kforvals
                        for r1tis=r1vals;
                            for flip = flipvals;
                                for Disp = Dispvals
                                    
                                    parms.mtis0(n) = mtis0;
                                    parms.f(n) = f;
                                    parms.cbva(n) = cbva;
                                    parms.bat(n) = bat;
                                    parms.bat2(n) = bat2;
                                    parms.kfor(n) = kfor;
                                    parms.r1tis(n) = r1tis;
                                    parms.flip(n) = flip;
                                    parms.Disp(n) = Disp;
                                    
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

%%  Testing by looking at the 10th entry 
n=10;
tp.mtis0 =parms.mtis0(n);
tp.f =    parms.f(n) ;
tp.cbva =  parms.cbva(n) ;
tp.bat =   parms.bat(n) ;
tp.bat2 = parms.bat2(n) ;
tp.kfor = parms.kfor(n);
tp.r1tis = parms.r1tis(n);
tp.flip = parms.flip(n) ;
tp.Disp = parms.Disp(n);

entry = gen_signals_160426(tp , timing_parms, 1, doSub);
%%

Nobs = length(entry);

dict = single(zeros(Ncombinations, Nobs, 'single'));
tic
%matlabpool open 4
parfor n=1:Ncombinations;
    dofigs = 0;
    if ~mod(n,1500)
        str = sprintf('echo Generating entry .... %d    of  %d   ',n, Ncombinations);
        system(str);
        dofigs = 1;
    end
    
    p = tp;
    
    p.mtis0 = parms.mtis0(n);
    p.f =    parms.f(n) ;
    p.cbva =  parms.cbva(n) ;
    p.bat =   parms.bat(n) ;
    p.bat2 = parms.bat2(n) ;
    p.kfor = parms.kfor(n);
    p.r1tis = parms.r1tis(n);
    p.flip = parms.flip(n) ;
    p.Disp = parms.Disp(n);

    entry = gen_signals_160426(p , timing_parms, dofigs, doSub);

    
    % Normalize all the entries of the dictionary
    tmp = entry - mean(entry);
    tmp = tmp / norm(tmp);
    entry = tmp;
    %
    dict(n,:) = entry;
    
end
%matlabpool close
toc
if dofigs
    figure(2);
    imagesc(real(dict) );  title('dictionary entries (rows)')
    drawnow
end

%save dictionary dict parms ;

end
