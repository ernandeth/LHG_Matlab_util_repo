function imgAscendFun(h, data, xyz, DesMat)
%function imgAscendFun(header, timeseries_data, [xyz_coords], DesMat)
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
% hard coded: choose whether we filter or pre-whiten

doWhiten = 0;
doFilter=1;

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

fprintf('\nBegin conditional probability Analysis ...');
ind = sub2ind([h.xdim, h.ydim, h.zdim], x,y,z);
exemplar = mean(data(:,ind),2);

if ~isempty(DesMat)
	% remove the known effects of interest (specified in design matrix):
	c = ones(1,size(DesMat,2));
	[t, beta_est, vBeta, RSS] = my_glm(DesMat,exemplar,c);
	exemplar = exemplar - DesMat*beta_est;
    fprintf('\nregressed out the design matrix')
end

MI = zeros(1,size(data,2));
asc = zeros(1,size(data,2));
H_asym = zeros(1,size(data,2));
rho = zeros(1,size(data,2));

if doWhiten
    % whitening the noise in the data
    exemplar = whiten(exemplar , DesMat);
end

% mask out pixels below this value
Thres1 =  max(data(10,:)) *0.1
warning off
for pix=1:size(data,2)
    if data(10,pix) > Thres1

        tseries = data(:,pix);
        if doWhiten
            tseries = whiten(tseries, DesMat);
        end
        
        if ~isempty(DesMat)
            c = ones(1,size(DesMat,2));
            [t, beta_est, vBeta,RSS] = my_glm(DesMat,tseries,c);
            tseries = tseries -  DesMat*beta_est;
        end
        
        if mod(pix,500)==0
            fprintf('\r Ascendency tests pix: %d of %d  ', pix, size(data,2))
        end

        [a, mi, H2, r] = causal03(exemplar',tseries',doFilter);

        asc(pix) = a;
        H_asym(pix) = H2;
        rho(pix)=r;
        MI(pix) = mi;
    end
end
warning on 

figure
% show histograms of the results.
subplot(221)
hist(asc(:),100);title('Ascendency')
%axis([-0.2 0.2 0 50])
subplot(222)
hist(H_asym(:),100);title('Histogram Asymmetry')
%axis([-0.2 0.2 0 50])
subplot(223)
hist(rho(:),100); title('Correlation COeff')
%axis([-1 1 0 50])
subplot(224)
hist(MI(:),100); title('Mutual Information')
%axis([-1 1 0 50])


asc = reshape(asc,h.xdim, h.ydim,h.zdim);
MI = reshape(MI,h.xdim, h.ydim,h.zdim);
H_asym = reshape(H_asym,h.xdim, h.ydim,h.zdim);
rho = reshape(rho,h.xdim, h.ydim,h.zdim);

global SPM_scale_factor
SPM_scale_factor = 1/1000;

write_hdr('MI.hdr',h);
write_hdr('Rho.hdr',h);
write_hdr('Asc.hdr',h);
write_hdr('H_asym.hdr',h);
write_img('MI.img',MI*1000,h);
write_img('Rho.img',rho*1000,h);
write_img('Asc.img',asc*1000,h);
write_img('H_asym.img',H_asym*1000,h);
save ascendency

return

