% function result = signal_lsq_func(estimate,paramters,data)
function result = signalDC_lsq(estimate,parms,data)

global diagn count_iter store % Define global variables
% disp('func called');

% Memory allocation
Bfield_est = estimate ; % Estimate
ind = parms.index;
% Bfield_est = estimate .* parms.weight;
% % With initial phase as a free parameter
% Bfield_amp =estimate(1:end-1);
% parms.initial_phase = estimate(end);
% [parms.phase_vec] = phase_calc(parms);



signal_estimate...
    = sum( (parms.m(ind) * ones(1,length(parms.phase_vec)) ...
    .* exp(-1i .* 2 .* pi .* (parms.yy_mat(ind) * parms.kro + parms.xx_mat(ind) * parms.kpe))...
    .* exp(-1i .* 2 .* pi .* parms.gambar .* Bfield_est * (cos(parms.sin_phase) .* ones(1,length(parms.phase_vec))) )), 1 );

if (nargin == 2)
    result = signal_estimate;
end

if (nargin == 3)
    %     signal_AC = data;
    %     RegTmp= zeros(parms.nv,parms.np/2);
    
    %     RegTmp(parms.index) = Bfield_est;
    %     [gx,gy] = gradient(abs(RegTmp),parms.xres,parms.yres);
    %     [gx , gy] = second_order_diff(abs(RegTmp) , parms);
    
    %     result =abs(signal_AC - signal_estimate).^2; % Objective function lsqnonlin
    result =(norm(abs(data -  signal_estimate))); %+ (parms.beta * norm(gx(:)+gy(:)));
    
    
    count_iter=count_iter+1;
    
    
    % Save important results
    if rem(count_iter,50) == 0
        
        diagn(store).Bf = Bfield_est;
        diagn(store).res = sum(result);
        diagn(store).sig2 = signal_estimate;
        
        
        store = store + 1;
        
    end
    
    
end

end


