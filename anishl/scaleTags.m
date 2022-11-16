function tags = scaleTags(min_c,npts_pr,intv,rangelen,comb,Nframes,tau,scan_time,tag_min)

global schedflag

coeffs = (min_c*ones(1,npts_pr))+ (intv* (my_ind2sub(rangelen*ones(1,npts_pr),comb)-1) ); % coeffs



% % for polynomial scheme
% 
% pol = abs(polyval(coeffs,1:Nframes)); %the polynomials
% init_tag = nan2zero(0.003*floor(fliplr( pol/(0.25*max(pol)) )/0.003)); % unscaled tags
% schedflag = 'Poly';
% % 



% for linInterp Scheme

t_series = zeros(1,Nframes);
spc = floor(Nframes/(npts_pr-1));
% t_series(1) = coeffs(1); t_series(end) = coeffs(end);
sample_points = [1 (1:(npts_pr-2))*spc Nframes];
t_series(sample_points) = coeffs;
init_tag = interp1(sample_points',t_series(sample_points),(1:Nframes)'); % Linear interpolation query
schedflag = 'LinInterp';



% total time and rescaling for scan_time
if (scan_time ~= Inf)
    tot_time = sum(init_tag);
    tags = ((scan_time-(tau*Nframes))*init_tag/tot_time) + tag_min; % scaled tags
    tags = nan2zero(tags);
else
    tags = init_tag;
end

tags = tags(:);