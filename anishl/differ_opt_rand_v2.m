function [diffSig] = differ_opt_rand_v2(params,delta,parm,tags,timing_parms,doSub)

%%


%%

parms.mtis0 =     1 ;
parms.Disp =      40;
parms.r1blood = 1/1.7;
parms.r1tis =     1/1.4  ;
parms.flip =      90*pi/180 ;


parms2.mtis0 =    1 ;
parms2.Disp =     40;
parms2.r1blood = 1/1.7;
parms2.r1tis =     1/1.4  ;
parms2.flip =      90*pi/180 ;
%%

cbva =      params(1) ;
bat_tis =   params(2);
bat_art =   params(3);
f=          params(4);
eta=        params(5) ;

frac = 0;

f = f*(1+frac);
cbva = cbva*(1+frac);
eta = eta*(1+frac);
bat_tis = bat_tis*(1+frac);
bat_art = bat_art*(1+frac);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +2.5% signal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_matrix = [cbva bat_tis bat_art f eta ];
p_matrix(parm) = (1+delta)*p_matrix(parm);

parms.cbva =      p_matrix(1); 
parms.f =         p_matrix(4); 

%%% PS MODEL %%%%
% parms.eta =       p_matrix(5);
% parms.bat_tis =   p_matrix(2); 
% parms.bat_art =   p_matrix(3);
%%%%%%%%%%%%%%%%%

%%% LHG MODEL %%%%
parms.kfor =  p_matrix(5);
parms.bat2 =  p_matrix(2);  
parms.bat =   p_matrix(3); 
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% -2.5% signal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_matrix = [cbva bat_tis bat_art f eta ];
p_matrix(parm) = (1-delta)*p_matrix(parm);

parms2.cbva =      p_matrix(1); 
parms2.f =         p_matrix(4);
%%% PS MODEL %%%%
% parms2.eta =       p_matrix(5);
% parms2.bat_tis =   p_matrix(2); 
% parms2.bat_art =   p_matrix(3);
%%%%%%%%%%%%%%%%%

%%% LHG MODEL %%%%
parms2.kfor =       p_matrix(5);
parms2.bat2 =   p_matrix(2);  
parms2.bat =   p_matrix(3); 
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig1 = (gen_signals_160426(parms,tags, timing_parms, 0,doSub));
sig2 = (gen_signals_160426(parms2,tags, timing_parms, 0,doSub));


p_matrix = [cbva bat_tis bat_art f eta];
diffSig = (sig1-sig2)/(2*delta*p_matrix(parm));