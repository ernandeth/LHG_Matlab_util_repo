%

%% optimize and generate 50 point ASL signal dictionary. The signal dictionary is optimized upon minimizing the correlation between each dictionary members of subtraction pairs.
% in order to minimize the impact of T1 and maximize the flow values, the
% dictionary generating function gen_dictionary_140712_r has fixed
% T1=1.4s.

% label_dur = x(1:20);
% PID = x(21:40);
% TR = x(41:60) ;
% using mean correlation+weight of total_sub+weight of tissue_sub
%Elapsed time is 4805.095320 seconds. G =0.6981
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;

ObjectiveFunction = @gmin_mix_0731_r; % traget function
nvars =50;    % Number of variables
% constrain the label durations:
random_l(1:nvars/2) = 0.9;
random_u(1:nvars/2) = 2.0;  
% constrain the PID s
random_l(nvars/2 + 1:nvars) = 0.01;
random_u(nvars/2 + 1:nvars) = 2.0;


LB = random_l;   % Lower bound
UB = random_u;  % Upper bound

tic

[x,fval,exitflag] = ga(ObjectiveFunction, nvars, [],[], [], [], LB,UB);


tvec = [x x ]; tvec(:) =0;
% duplicate the timing vector so that we can do teh control and the tag the
% same
tvec(1:2:end) = x;
tvec(2:2:end) = x;

%keyboard
timing_parms.label_dur = tvec(1:end/2)';
timing_parms.PID =       tvec(end/2 +1 :  end)';
timing_parms.TR =        timing_parms.PID + timing_parms.label_dur + 0.07;


% [dict, parms] = gen_dictionary_140724(timing_parms,0);
% entry = gen_signals_140707(parms(300) , timing_parms, 1,0);

[dict,diff_tissue, parms] = gen_dictionary_140712_r(timing_parms,1);

save optimalTiming_140731c timing_parms dict parms

toc

figure(2)
subplot(211)
imagesc(dict)
title('dictionary members');colorbar;

xc = (corrcoef(dict'));
subplot(212)
imagesc(abs(xc))
title('correlation map');colorbar;

xc = abs(xc(:));
G = mean(xc)

quit

% data=dict-repmat(mean(dict,1),size(dict,1),1);
% [W, EvalueMatrix] = eig(cov(data'));
% Evalues = diag(EvalueMatrix);
% Evalues = Evalues(end:-1:1);
% W = W(:,end:-1:1); W=W';  


% %%% investigate flow susceptibility of timing_parms

% testdict=zeros(100,nvars);
% testfvals=linspace(0,100,100)/6000;
% testmtis0 =     1 ;
% testcbva =      0.02 ;
% testtransit =   1.2 ;
% testkfor =      0.2; % 1e-2 ;
% testr1tis =     1/1.4  ;
% testbeta =      90*pi/180 ; % flip angle in radians
% testDisp =      30;
% figure(9)
% for n=1:length(testfvals)
%     testparms(n).f=testfvals(n);
%     testparms(n).mtis0 = testmtis0;
%     testparms(n).cbva = testcbva;
%     testparms(n).transit = testtransit;
%     testparms(n).kfor = testkfor;
%     testparms(n).r1tis = testr1tis;
%     testparms(n).beta = testbeta;
%     testparms(n).Disp = testDisp;
%     testdict(n,:)=gen_signals_140707(testparms(n),timing_parms,0,0);
%     n=n+1;
% end
% subplot(211)
% imagesc(testdict);colorbar;
% title('generated signal upon different CBF values')
% subplot(212)
% testxc=corrcoef(testdict');
% imagesc(testxc);colorbar;
% title('Pierson Correlation Coefficient Matrix')

