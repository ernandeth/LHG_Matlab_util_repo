% a script to look at the sensitivity of different timing schedules to each
% parameter.  

%Use partial derivatives and the cramer rao bound



parms = struct( ...
    'mtis0', 1,...
    'f', 0.01 , ...
    'cbva' ,  0.01 , ...
    'bat', 1.0,...
    'bat2', 1.2,...
    'kfor', 0.02, ...
    'r1tis', 0.8,...
    'flip', deg2rad(75), ...
    'Disp', 0); % 40

close all
!rm sensitivities.txt CRLBs.txt
pfs = fopen('sensitivities.txt' , 'wt') 
pfc = fopen('CRLBs.txt' , 'wt') 
%
% creating ASL timing schedule that could be used to fit transit time
Nframes = 50;

tmp = linspace(0.08, 3, Nframes/2);
t_tag = zeros(Nframes,1 );
t_tag(1:2:end) = tmp;
t_tag(2:2:end) = tmp;

t_delay =  ones(size(t_tag)) * 0.08;
isLabel = ones(size(t_tag));
isLabel(1:2:end) = 0;
t_aq = ones(size(t_tag)) * 0.035;
t_adjust = 5 - t_tag - t_delay - t_aq;

timing_parms.Nlabel_group = 1;
timing_parms.t_tag = t_tag;
timing_parms.t_delay = t_delay;
timing_parms.t_adjust = t_adjust;
timing_parms.isLabel = isLabel;
timing_parms.order = 1;
timing_parms.t_aq = t_aq;

fp_sensitivity_analysis
set(gcf,'Name', 'Linear increase')

 fprintf(pfc, '%.3g\t',100* sqrt(diag(CRLB))' ./ p' );
 fprintf(pfc, '\n' );
 fprintf(pfs, '%.3g\t', sensitivity');
 fprintf(pfs, '\n' );
 

%%  Now make a schedule for the perfusion part of the ASL experiment

Nframes = 50;


t_tag = ones(Nframes,1 ) * 1.8;
t_delay =  ones(size(t_tag)) * 1.8;
isLabel = ones(size(t_tag));
isLabel(1:2:end) = 0;
isLabel(1:8) = 0;
t_aq = ones(size(t_tag)) * 0.035;
t_adjust = ones(size(t_tag)) * 4  ;
t_adjust = t_adjust - t_tag - t_delay - t_aq;

timing_parms.Nlabel_group = 1;
timing_parms.t_tag = t_tag;
timing_parms.t_delay = t_delay;
timing_parms.t_adjust = t_adjust;
timing_parms.isLabel = isLabel;
timing_parms.order = 1;
timing_parms.t_aq = t_aq;

fp_sensitivity_analysis
set(gcf,'Name', 'standard PCASL')

 fprintf(pfc, '%.3g\t',100* sqrt(diag(CRLB))' ./ p' );
 fprintf(pfc, '\n' );
 fprintf(pfs, '%.3g\t', sensitivity');
 fprintf(pfs, '\n' );
%
%% other schedules
Npoints = 500;
timing_parms.t_delay =  0.05 * ones(1,Npoints) ;
timing_parms.t_adjust = 0.05 * ones(1,Npoints) ;
timing_parms.isLabel = round(rand(1,Npoints))  ;
timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
timing_parms.t_aq = ones(1,Npoints) * 0.035;

supercontrols = randperm(Npoints);
supercontrols=supercontrols(1:round(Npoints/4));
timing_parms.isLabel(supercontrols) = -1;

timing_parms.t_tag = 1*abs(sinc(linspace(-2, 1, Npoints))) + 0.05 ;
fp_sensitivity_analysis
set(gcf,'Name', 'sinc')

 fprintf(pfc, '%.3g\t',100* sqrt(diag(CRLB))' ./ p' );
 fprintf(pfc, '\n' );
 fprintf(pfs, '%.3g\t', sensitivity');
 fprintf(pfs, '\n' );
 
timing_parms.t_tag = 3*abs(sinc(linspace(-2, 1, Npoints))) + 0.05 ;
fp_sensitivity_analysis
set(gcf,'Name', '3 x sinc')

 fprintf(pfc, '%.3g\t',100* sqrt(diag(CRLB))' ./ p' );
 fprintf(pfc, '\n' );
 fprintf(pfs, '%.3g\t', sensitivity');
 fprintf(pfs, '\n' );

timing_parms.t_tag = 2*abs((linspace(0, 1, Npoints) .^2)) + 0.05 ;
fp_sensitivity_analysis
set(gcf,'Name', 'quadratic')

 fprintf(pfc, '%.3g\t',100* sqrt(diag(CRLB))' ./ p' );
 fprintf(pfc, '\n' );
 fprintf(pfs, '%.3g\t', sensitivity');
 fprintf(pfs, '\n' );


timing_parms.t_tag = 1-1*abs((linspace( 1,0, Npoints) .^2)) + 0.05 ;
fp_sensitivity_analysis
set(gcf,'Name', 'reverse quadratic')

 fprintf(pfc, '%.3g\t',100* sqrt(diag(CRLB))' ./ p' );
 fprintf(pfc, '\n' );
 fprintf(pfs, '%.3g\t', sensitivity');
 fprintf(pfs, '\n' );
 
timing_parms.t_tag = 2*abs((linspace(0, 1, Npoints) .^3)) + 0.05 ;
fp_sensitivity_analysis
set(gcf,'Name', 'cubic')

 fprintf(pfc, '%.3g\t',100* sqrt(diag(CRLB))' ./ p' );
 fprintf(pfc, '\n' );
 fprintf(pfs, '%.3g\t', sensitivity');
 fprintf(pfs, '\n' );
 
timing_parms.t_tag = [0.5*triang(Npoints/2) ;  1*triang(Npoints)]' + 0.05 ;   
timing_parms.t_tag = timing_parms.t_tag(1:Npoints);
fp_sensitivity_analysis
set(gcf,'Name', 'triangle')

 fprintf(pfc, '%.3g\t',100* sqrt(diag(CRLB))' ./ p' );
 fprintf(pfc, '\n' );
 fprintf(pfs, '%.3g\t', sensitivity');
 fprintf(pfs, '\n' );
 
% timing_parms.t_tag = 2*abs(sin(linspace(0,4*pi, Npoints))) + 0.05 ;
% fp_sensitivity_analysis
% set(gcf,'Name', 'sin')
% 
%  fprintf(pfc, '%.3g\t',100* sqrt(diag(CRLB))' ./ p' );
%  fprintf(pfc, '\n' );
%  fprintf(pfs, '%.3g\t', sensitivity');
%  fprintf(pfs, '\n' );
%  
fclose(pfs);
fclose(pfc);

%%
tmp = 2.5*abs(sinc(linspace(-5, 0 , Npoints) ));
timing_parms.t_tag = tmp + tmp(end:-1:1) +0.05;

fp_sensitivity_analysis
set(gcf,'Name', '2.5xsinc .. 5lobe')

sum(timing_parms.t_tag)
