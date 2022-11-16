Npoints = 300;
durs = 2*abs( (linspace(1, 0.1,Npoints)) .* cos(linspace(0,4*pi,Npoints))) + 0.02 ;
durs = 2*abs( sinc(linspace(-2,2,Npoints))) + 0.02 ;
timing_parms.t_tag =  durs ;
timing_parms.t_delay = 0.05 * ones(1,Npoints) ;
timing_parms.t_adjust = 0.05 * ones(1,Npoints) ;
timing_parms.isLabel = round(1.5-3*rand(1,Npoints))  ;
timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
timing_parms.t_aq = ones(1,Npoints) * 0.035;
 

%% synth. training parameters (targets)
parms0 = struct( ...
    'mtis0', 1,...
    'f', 0.01 , ...
    'cbva' ,  0.01 , ...
    'bat', 1.0,...
    'bat2', 1.2,...
    'kfor', 0.5, ...
    'r1tis', 0.8,...
    'flip', pi/4, ...
    'Disp', 50);

global dt
dofigs = 0
doSub=0
tmp =[];
for dt=[0.1:0.2:2]*1e-3

    entry = gen_signals_160426(parms0 , timing_parms, dofigs, doSub);
    tmp = [tmp entry];
end

delts = repmat(tmp(:,1), 1,size(tmp,2)) -tmp;
for n=1:size(tmp,2)
    plot(delts(:,n))
    hold on
  
end
