function G = gmin_mix_0731_r(x)

% used to adjust weighting parameters of ma and ta
% this function  takes in a vector of timing parameters and uses it to
% create a dictuinary of ASL signals.

% in order to minimize the impact of T1 and maximize the flow values, the
% dictionary generating function gen_dictionary_140712_r has fixed
% T1=1.4s.

%
% label_dur = x(1:20);
% PID = x(21:40);
% TR = x(41:60) ;
%
% it returns the MSE of the correlation matrix of the dictionary
%
gamma=.005;
lambda=.005;
% phi=2;
% nvars=50;

tvec = [x x]; tvec(:) =0;
% duplicate the timing vector so that we can do teh control and the tag the
% same
tvec(1:2:end) = x;
tvec(2:2:end) = x;


%keyboard
timing_parms.label_dur = tvec(1:end/2)';
timing_parms.PID =       tvec(end/2 +1 :  end)';
timing_parms.TR =        (timing_parms.PID + timing_parms.label_dur + 0.070);
% TR includes just enough time for a single slice and 35 ms. after
% acquisition
% 
% [dict, parms] = gen_dictionary_140711(timing_parms,0);
[diff_tot,diff_tissue, parms] = gen_dictionary_140712_r(timing_parms,1);

%entry = gen_signals_140420(parms(1) , timing_parms, 1);
% matlabpool('open',4);
% diff_tot=mean(diff_tot(:));
% diff_tissue=mean(diff_tissue(:));
% matlabpool('close');


% % PCA analysis
% data=dict-repmat(mean(dict,1),size(dict,1),1);
% [W, EvalueMatrix] = eig(cov(data'));
% Evalues = diag(EvalueMatrix);
% Evalues = Evalues(end:-1:1);
% W = W(:,end:-1:1); W=W';  
% dim=find(Evalues<0.001,1);
% Sreal=std(Evalues(2:dim-1));


%correlation computation
xc = (corrcoef(diff_tot'));
xc = abs(xc(:));
Greal=mean(xc)


diff_tot=mean(diff_tot(:));
diff_tissue=mean(diff_tissue(:));

% G=-Sreal-gamma*10*log(diff)/log(10)-lambda*log(dim)/log(10)
G=Greal-gamma*10*log(diff_tot)/log(10)-lambda*10*log(diff_tissue)/log(10);

% %% rewarding variations in cbf flow
% testdict=zeros(100,nvars);
% testfvals=linspace(0,100,100)/6000;
% testmtis0 =     1 ;
% testcbva =      0.02 ;
% testtransit =   1.2 ;
% testkfor =      0.2; % 1e-2 ;
% testr1tis =     1/1.4  ;
% testbeta =      90*pi/180 ; % flip angle in radians
% testDisp =      30;
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
% testxc=corrcoef(testdict');
% 
% Fscore=std(testxc(:));

% G=Greal-gamma*10*log(diff)/log(10)+phi*(1-Fscore)


return
    
    
