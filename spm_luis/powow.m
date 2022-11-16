% function N = power_ranger_jr(x, alpha, bf, error, power_threshold)
% function N = power_ranger_jr(x, alpha, bf, error, power_threshold)
%
% computes the Number of measurements
% to achieve a given power for a t-test accross
% several images of every voxel.
%
% Luis Hernandez
% University of Michigan
% Last Edit 2 - 15 -2000
%
% The null hypothesis is that each images has a value equal
% to the mean value of all the images in the group.
% The alternative hypothesis is that they are not (double sided test)
% a p-value less than alpha means that they are different. 
%
% This function computes the probability of being right if we reject the 
% null hypothesis of being equal, based on a given effect size.
%
% Arguments and defaults:
%   x:       data array
% 	alpha: 	significance level (0.05)
% 	bf :   	Bonferroni correction factor. 
%   error:   Effect size to be assumed, when calculating the power (5%)
%   power_threshold:       Desired power level
%
% 


  %x = [130 131 131 120 133 125 135 131 131 131 132 132 140 150 100]
  bf = 1    
  alpha = 0.05
  error = sqrt(var(x))*1
  power_threshold = 0.8;
  error = 10^5-10^3;

alpha = alpha/(bf*2);  % 2-tailed Significance level after Bonferroni correction
N = 3; %x_size;
power = 0;

while power < power_threshold



    % compute statistics for all images
    x_mean  = 10^5;
    x_sigma = 0.1*10^5/N;
    %(var(data))'/(N_correction);
    x_var = x_sigma^2;
    
    

    % make a plot of the normal ditribution
%    f = normpdf(x,x_mean, x_sigma);
%    hold off
%    plot(x,f,'b*')
%    hold on
%    xx = [x_mean-2*x_var:1:x_mean + 2*x_var];
%    ff = normpdf(xx,x_mean, x_sigma);
%    plot(xx,ff,'b');
%    %
%
%
%
%    % compute ditribution assuming that we are wrong and the truth is one 
%    % some standard deviations higher. The t_score in this ideal distribution
%    % is then given by
%
%    f = normpdf(x, x_mean + error , x_sigma);
%    plot(x, f,'g*')
%
%    x_ideal = x + error ;
%    xi = [min(x_ideal):0.1:max(x_ideal)];
%    fi = normpdf(xi, x_mean +  error , x_sigma);
%    plot(xi,fi,'g')
%
%    %
	
	
    % find critical score for this significance level(alpha)
    t = 0;
    p=100;  %(some initial value)
    while (p >= alpha),
%  		p=(1-normcdf(t,x_mean, x_sigma));
  		p=(1-spm_Ncdf(t,x_mean, x_sigma));
  		t = t + x_sigma/10;
      end
      x_alpha = t;
      alpha;

      %plot(x_alpha, normpdf(x_alpha, x_mean, x_sigma), 'r*')
      %plot(x_alpha, normpdf(x_alpha, x_mean + error, x_sigma), 'r*')
      plot(x_alpha, spm_Npdf(x_alpha, x_mean, x_sigma), 'r*')
      plot(x_alpha, spm_Npdf(x_alpha, x_mean + error, x_sigma), 'r*')
      
      axis([min(x)-10 max(x)+10 + error 0 0.1]);
      drawnow

      
      beta =  spm_Ncdf(x_alpha, x_mean + error, x_sigma);
      power = 1 -beta
      x_sigma
            
      if power < power_threshold
          N = N+1
      end
      

  end
  
  
return
