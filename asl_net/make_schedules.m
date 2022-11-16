
addpath  /home/hernan/matlab/flow/fingerprint

Npoints = 200
Npoints = 400

timing_parms.t_delay =  0.05 * ones(1,Npoints) ;
timing_parms.t_adjust = 0.05 * ones(1,Npoints) ;
timing_parms.isLabel = round(rand(1,Npoints))  ;
supercontrols = randperm(Npoints);
supercontrols= supercontrols(1:round(Npoints/3));
timing_parms.isLabel(supercontrols) = -1;

timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
timing_parms.t_aq = ones(1,Npoints) * 0.035;

timing_parms.t_tag = 3*abs(sinc(linspace(2, 0, Npoints))) + 0.05 ;
timing_parms.t_tag = 2*abs((linspace(0, 1, Npoints) .^2)) + 0.05 ;
timing_parms.t_tag = 1-1*abs((linspace( 1,0, Npoints) .^2)) + 0.05 ;
timing_parms.t_tag = 1*abs(sinc(linspace(-2, 0, Npoints))) + 0.05 ;
%timing_parms.t_tag = [0.5*triang(Npoints/2) ;  triang(Npoints)]' + 0.05 ;   timing_parms.t_tag = timing_parms.t_tag(1:Npoints);
timing_parms.t_tag = 1.8*abs(sin(linspace(0,4*pi, Npoints))) + 0.05 ;

totalsecs = sum( timing_parms.t_tag) + sum( timing_parms.t_delay) + sum( timing_parms.t_adjust) + sum( timing_parms.t_aq)

parms0 = struct( ...
    'mtis0', 1,...
    'f', 0.008 , ...
    'cbva' ,  0.01 , ...
    'bat', 1.0,...
    'bat2', 1.2,...
    'kfor', 0.015, ...
    'r1tis', 0.8,...
    'flip', pi/2, ...
    'Disp', 20);


% show the signal for one case
doSub = 0;
dofigs = 1;

figure


entry = gen_signals_180320(parms0, timing_parms, dofigs, doSub);
var(entry)
    
parms_dummy = struct( ...
    'mtis0', 1,...
    'f', 1 , ...
    'cbva' , 1 , ...
    'bat', 1,...
    'bat2', 1,...
    'kfor', 1, ...
    'r1tis', 1,...
    'flip', 1, ...
    'Disp', 1);
ders = calc_partials(timing_parms, parms0) ;


figure
subplot(211)
plot(timing_parms.t_tag)
ylabel('label duration')
xlabel('Scan')
fatlines

subplot(212)
plot(entry)
ylabel('Signal')
xlabel('Scan')
fatlines

print -dpng ASL_MRF_signal

t_tag = timing_parms.t_tag';
t_delay = timing_parms.t_delay';
t_adjust = timing_parms.t_adjust';
isLabel = timing_parms.isLabel';

!mkdir timing_sinusoid3_400

save timing_sinusoid3_400/t_tags.txt t_tag -ascii
save timing_sinusoid3_400/t_adjusts.txt t_adjust -ascii
save timing_sinusoid3_400/t_delays.txt t_delay -ascii
save timing_sinusoid3_400/labelcontrol.txt isLabel -ascii

%%
parms = parms0;
test_vals = 0.01: 0.05:0.5;
diffs = [];
     
t_delay = 0.5* rand(1,400) + 0.05;
 tmp =  round(t_delay * 1000/20) * 20/1000;
 timing_parms.t_delay  = tmp;
     
for n=1:length(test_vals)
     %t_tag= 2 * abs(sinc(linspace(-1, 2 , Npoints))) + 0.05 + test_vals(n) ;
     t_tag= 1.8 * abs(sin(linspace(-pi, pi , Npoints))) + 0.05 + test_vals(n) ;
     %t_tag = 2 * (abs(linspace( 0, 1, Npoints) .^ 2 )) + 0.05 + test_vals(n);
     % t_tag =2.5 *  [triang(Npoints/2) ;  0.3*triang(Npoints/2)]' + 0.05 ; 
     % t_tag = timing_parms.t_tag(1:Npoints);

%      load /home/hernan/matlab/anishl/minLinInterp600LH2_700mod4_20_TR_fix.mat
%      timing_parms = timing_parms_min;
%      t_tag = timing_parms.t_tag;
     
     % round off to 20 ms. intervals
     tmp =  round(t_tag * 1000/20) * 20/1000;
     %timing_parms.t_tag  = [tmp(1:2:end)   tmp(2:2:end)  ];
     timing_parms.t_tag  = tmp;
     

    
     %timing_parms.t_delay(:)  = 0.05;

     totalsecs = sum( timing_parms.t_tag) + sum( timing_parms.t_delay) + sum( timing_parms.t_adjust) + sum( timing_parms.t_aq)
     entry = gen_signals_180320(parms, timing_parms, 1, 0);

    drawnow
    
    % Look at the sensitivity to phys parms:
    % calculate dSignal/dParm
    ders = calc_partials(timing_parms, parms) ;
    
    % order of derivatives from the function:
    diffs(n).mtis0   = ders(1,1) ;
    diffs(n).f           = ders(2,2);
    diffs(n).cbva     = ders(3,3);
    diffs(n).bat       = ders(4,4);
    diffs(n).bat2     = ders(5, 5) ;
    diffs(n).kfor      = ders(6,6);
    diffs(n).r1tis     = ders(7,7);
    diffs(n).flip       = ders(8,8);
    diffs(n).Disp     = ders(9,9);
    
end

diffs

dS_df = zeros(1, length(test_vals));
dS_dbat = zeros(1, length(test_vals));
dS_dcbva = zeros(1, length(test_vals));
dS_dkfor = zeros(1, length(test_vals));


for n=1:length(test_vals),     
    dS_df(n) = ( diffs(n).f); 
    dS_dbat(n) = ( diffs(n).bat); 
    dS_dcbva(n) = ( diffs(n).cbva);
    dS_dkfor(n) = ( diffs(n).kfor);
end

subplot(221), plot(dS_df), title('dS/df')
subplot(222), plot(dS_dbat), title('dS/dbat')
subplot(223), plot(dS_dcbva), title('dS/dcbva')
subplot(224), plot(dS_dkfor), title('dS/dkfor')

% 
% load /home/hernan/matlab/anishl/minLinInterp600LH2_700mod4_20_TR_fix.mat
% timing_parms = timing_parms_min;
% ders = calc_partials(timing_parms, parms) ;
% gen_signals_180320(parms, timing_parms, 1, 0);
% 


t_tag = timing_parms.t_tag';
t_delay = timing_parms.t_delay';
t_adjust = timing_parms.t_adjust';
isLabel = timing_parms.isLabel';

!mkdir timing_sinusoid2del_400

save timing_sinusoid2del_400/t_tags.txt t_tag -ascii
save timing_sinusoid2del_400/t_adjusts.txt t_adjust -ascii
save timing_sinusoid2del_400/t_delays.txt t_delay -ascii
save timing_sinusoid2del_400/labelcontrol.txt isLabel -ascii
