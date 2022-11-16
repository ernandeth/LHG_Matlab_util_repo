function kern = T2filter(type,raw,scaninfo,doSHOWPIX,fh)
% T2FILTER: function for detrending echo train based on T2 decay
%|
%|  Accepts "b-matrix" (raw data reshaped into {# data pts per echo} x
%|      {# of shots} x {# of slices}) and returns version of b corrected
%|      for T2 decay and fitted T2 (units of samples) and x0 for purpose of
%|      plotting
%|
%|  Parameters:
%|      type: specify how T2 values act on data
%|          1: Each echo train (frame/coil/leaf) has an individual T2 fit
%|          2: Each frame/coil has an individual T2 fit
%|          3: One T2 fit for the entire dataset
%|      raw: raw data
%|      scaninfo: scaninfo struct from read_raw_3d
%|      doSHOWPIX: specify if you would like to show the fit figure
%|      fh: figure handle if doSHOWPIX is on
%|
    
    % Get info from scaninfo struct:
    TE = scaninfo.te;
    ts = scaninfo.ts;
    nslices = scaninfo.nslices;
    nleaves = scaninfo.npr;
    ndat = scaninfo.ndat;
    ncenter = scaninfo.ncenter;
    nframes = scaninfo.nphases;
    ncoils = scaninfo.ncoils;

    ndat_all = round(TE*1e-3/ts);
    navs = round(ndat/2-ncenter/2):round(ndat/2+ncenter/2);

    % Set defaults:
    if nargin<4 || isempty(doSHOWPIX)
        doSHOWPIX = 1;
    end
    
    if (nargin<5 || isempty(fh)) && doSHOWPIX
        fh = figure;
    elseif doSHOWPIX
        figure(fh);
    end

    % Initialize important arrays:
    R2s = zeros(nleaves,nframes,ncoils);
    x = reshape(1:ndat_all*nslices,ndat_all,nslices);
    kern = zeros([nframes ndat nleaves nslices ncoils]);

    for iframe = 1:nframes
        for icoil = 1:ncoils
            b = reshape(raw(iframe,:,icoil),[ndat,nleaves,nslices]);
            b_all = cat(1,b,NaN(ndat_all-ndat,nleaves,nslices));
            for ileaf = 1:nleaves
                y_navs = mean(abs(b_all(navs,ileaf,:)),1);
                x_navs = x(round(ndat/2),:);

                A = x_navs(:).^[1 0];
                coeffs = (A'*A)^(-1) * A' * log(y_navs(:));
                R2s(ileaf,iframe,icoil) = coeffs(1);
                s0 = exp(coeffs(2));

                t2fit = s0*exp(R2s(ileaf,iframe,icoil)*(ndat_all*floor(x/ndat_all) + ndat/2));

                if doSHOWPIX
                    subplot(2,1,1)
                    plot(x,squeeze(abs(b_all(:,ileaf,:)))), hold on
                    plot(x(:),t2fit(:)), hold off
                    xlabel('Data points'), ylabel('Magnitude')
                    ylim([0 max(abs(raw),[],'all')])
                    title(sprintf('T2 fit for frame %d, coil %d, leaf %d\nEstimated T2: %.1fms',...
                        iframe,icoil,ileaf,abs(TE/(ndat*R2s(ileaf,iframe,icoil)))))
                    hold off
                    subplot(2,1,2)
                    histogram(abs(TE./(ndat*R2s(R2s~=0))),'BinWidth',25)
                    title(sprintf('Histogram of fitted T2 values\nMean: %.1fms STD: %.1fms',...
                        mean(abs(TE./(ndat*R2s(R2s~=0))),'all'), std(abs(TE./(ndat*R2s(R2s~=0))))))

                    xlabel('T2 value (ms)'), ylabel('Counts')
                    drawnow
                end
            end
        end
    end

    switch type
        case 1
        case 2
            for iframe = 1:nframes
                for icoil = 1:ncoils
                    R2s(:,iframe,icoil) = mean(R2s(:,iframe,icoil),1);
                end
            end
        case 3
            R2s(:,:,:) = mean(R2s,'all');
        otherwise
            error('Invalid type: %s\n',string(type))
    end

    for iframe = 1:nframes
        for icoil = 1:ncoils
            for ileaf = 1:nleaves
                t2fit = exp(R2s(ileaf,iframe,icoil)*(ndat_all*floor(x/ndat_all) + ndat/2));
                kern(iframe,:,ileaf,:,icoil) = t2fit(1:ndat,:);
            end
        end
    end

end

