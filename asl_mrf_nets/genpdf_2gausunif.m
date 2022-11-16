function pdf = genpdf_2gausunif(mean1,mean2,sigma1,sigma2,support)

pt1 = gaussmf(support,[sigma1 mean1]); % Gaussian 1
pt2 = gaussmf(support,[sigma2 mean2]); % Gaussian 2
% pt3 = (support(end)/100.0) * (ones(size(support))*(support(end)-support(1))^-1); % Uniform 
pt3 = 0;

pdf = pt1 + pt2 + pt3; % sum over support
pdf = pdf/(sum(pdf)); % normalizing to sum 1