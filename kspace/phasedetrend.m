function b_corr = phasedetrend(b,Ncenter,PDtype,spec,doSHOWPIX,fh)
%PHASEDETREND: function for correcting phase by evaluating drift at
%   center of kspace
%|
%|  Accepts "b-matrix" (raw data reshaped into {# data pts per echo} x
%|      {# of shots} x {# of slices}) and returns version of b corrected
%|      for phase drift
%|
%|  Parameters:
%|      b: "b-matrix" from spiralrecon
%|      Ncenter: number of points at center of kspace
%|      PDtype: 'lsqNav' for poly detrend based on navigator points in echo
%|          'lsqAll' for poly detrend based on all points in echo
%|          'HPF' for high pass filtering
%|      spec: max order if using poly detrend, normalized cutoff frequency 
%|          if using hpf (expressed in frac of nyquist spectrum)
%|      doSHOWPIX: option to display plots for debugging
%|      fh: figure handle to use (default generates a new figure)
%|

    % Prep type comparison:
    PDtype = lower(PDtype);
    PDtype_int = 1*strcmp(PDtype,'lsqnav') + ...
        2*strcmp(PDtype,'lsqall') + ...
        3*strcmp(PDtype,'hpf');
    
    % Setting defaults:
    switch PDtype_int
        case 0 % Check that PDtype is valid
            error("Invalid PDtype: %s",PDtype)
        case 1 || 2
                if nargin<4 || isempty(spec), spec = 0; end
        case 3
            if nargin<4 || isempty(spec), spec = 20; end
    end
    
    if nargin<5 || isempty(doSHOWPIX)
        doSHOWPIX = 1;
    end
    
    if nargin<6 || isempty(fh)
        fh = figure;
    end
    
    % Determine navigating points:
    allpts = (1:size(b,1))';
    navpts = (round(size(b,1)/2 - Ncenter/2 + 1) : round(size(b,1)/2 + Ncenter/2))';
    
    % Initialize design matrix:
    n = spec:-1:0;
    switch PDtype_int
        case 1 % lsq based on navigator points
            A = navpts.^n;
        case 2 % lsq based on all points
            A = allpts.^n;
        case 3 % low pass filter (no lsq)
            A = [];
    end

    % Initialize corrected echo matrix:
    b_corr = zeros(size(b));
    b_phase = zeros(size(b));
    sub_phase = zeros(size(b));
    b_phase_corr = zeros(size(b));

    for shot = 1:size(b,2)
        for slice = 1:size(b,3)
            % Unwrap phase of echo:
            b_phase(:,shot,slice) = unwrap(angle(b(:,shot,slice)));
            
            if PDtype_int == 1 || PDtype_int == 2
                % Determine regressors using least squares:
                if PDtype_int == 1
                    coeffs_phase = (A'*A)^(-1) * A' * b_phase(navpts,shot,slice);
                else
                    coeffs_phase = (A'*A)^(-1) * A' * b_phase(allpts,shot,slice);
                end
                sub_phase(:,shot,slice) = (1:length(b(:,shot,slice)))'.^(spec:-1:0) * coeffs_phase;
                
                % Force correction to have same 0-crossing:
                sub_phase(:,shot,slice) = sub_phase(:,shot,slice) ...
                    - (sub_phase(round(end/2),shot,slice) - b_phase(round(end/2),shot,slice));
                
                % Calculate new echo phase with lsq correction:
                b_phase_corr(:,shot,slice) = b_phase(:,shot,slice) - sub_phase(:,shot,slice);
            else
                % Calculate highpass filtered data:
                fspec = linspace(-1/2,1/2,size(b,1))';
                KERNEL = 1/2*fftshift(tanh(200*pi*(fspec-spec/4)) ...
                    + tanh(200*pi*(-fspec-spec/4))) + 1;
                
                % Calculate new echo phase with highpass filter correction:
                b_phase_corr(:,shot,slice) = real(ifft(fft(b_phase(:,shot,slice)).*KERNEL));
            end
            
            % Calculate corrected echo:
            b_corr(:,shot,slice) = abs(b(:,shot,slice)) .* ...
                exp( sqrt(-1) * b_phase_corr(:,shot,slice) );
        end
    end

    
    % Create figure for debugging:
    if doSHOWPIX
        fspec = fh;
        fspec.Position(3:4) = fspec.Position(1:2) + [1400 500];
        shot = 1;
        for slice = 1:size(b,3)
            % Plot entire uncorrected vs. corrected echo:
            subplot(2,3,1)
                plot(allpts,b_phase(allpts,shot,slice),'-r')
                hold on
                plot(allpts,b_phase_corr(allpts,shot,slice),'-b')
                % If poly detrend, plot the fits as well:
                if PDtype_int == 1 || PDtype_int == 2
                    plot(allpts,sub_phase(allpts,shot,slice),'--r')
                    plot(allpts,0*allpts,'--b')
                end
                hold off
                
                % Properties:
                b_rev = [b_phase(allpts,shot,slice);
                    b_phase_corr(allpts,shot,slice)];
                dev = 1.1*max(abs(b_rev - mean(b_rev,'all')));
                bounds = mean(b_rev) + dev*[-1 1];
                ylim(bounds), yticks([bounds(1) bounds(end)]), ylabel('Radians')
                yticklabels([sprintf("%.2fπ",bounds(1)) sprintf("%.2fπ",bounds(2))])
                xlabel("All data points in echo")
                legend('Uncorrected phase', 'Corrected phase',...
                    'Location','NorthOutside','Orientation','horizontal')
                legend('boxoff'), ylabel('Radians')
                title('Phase detrending for first echo')

            % Plot uncorrected vs. corrected echo @ navigator points:
            subplot(2,3,4)
                % Uncorrected phase on left yaxis:
                    yyaxis left
                    plot(navpts,b_phase(navpts,shot,slice),'-r')
                    hold on
                    % If poly detrend, plot the fits as well:
                    if PDtype_int == 1 || PDtype_int == 2
                        plot(navpts,sub_phase(navpts,shot,slice),'--r')
                    end
                    hold off
                    
                    % Properties:
                    b_rev = b_phase(navpts,shot,slice);
                    dev = 1.1*max(abs(b_rev - mean(b_rev,'all')));
                    bounds = mean(b_phase(navpts,shot,slice)) + dev*[-1 1];
                    ylim(bounds), yticks([bounds(1) bounds(2)])
                    ylabel('Radians')
                    yticklabels([sprintf("%.2fπ",bounds(1)) sprintf("%.2fπ",bounds(2))])

                % Corrected phase on right yaxis:
                    yyaxis right
                    plot(navpts,b_phase_corr(navpts,shot,slice),'-b')
                    hold on
                    % If poly detrend, plot the fits as well:
                    if PDtype_int == 1 || PDtype_int == 2
                        plot(navpts,0*navpts,'--b')
                    end
                    hold off
                
                    % Properties:
                    bounds = mean(b_phase_corr(navpts,shot,slice)) + dev*[-1 1];
                    ylim(bounds), yticks([bounds(1) bounds(2)]), ylabel('Radians')
                    yticklabels([sprintf("%.2fπ",bounds(1)) sprintf("%.2fπ",bounds(2))])
                    xlabel("Navigating data points in echo")

                % Set axis colors:
                ax = gca;
                ax.YAxis(1).Color = 'r';
                ax.YAxis(2).Color = 'b';

            % Surf plot uncorrected vs. corrected echos of each slice:
            subplot(2,3,[2 3 5 6])
                surf(squeeze(b_phase(navpts,shot,:)), ...
                    'FaceAlpha', 0.5), hold on
                surf(squeeze(b_phase_corr(navpts,shot,:)), ...
                    'FaceAlpha', 0.5,'FaceColor','red')
                hold off

                % Properties:
                legend('Uncorrected phase', 'Corrected phase', ...
                    'Location','SouthOutside','Orientation','horizontal')
                legend('boxoff')
                xlabel('Slice index'), ylabel('Data points')
                yticks(linspace(1,length(navpts),5))
                yticklabels(linspace(navpts(1),navpts(end),5))
                b_rev = [b_phase(navpts,shot,:);
                    b_phase_corr(navpts,shot,:)];
                dev = 1.1*max(abs(b_rev - mean(b_rev,'all')),[],'all');
                bounds = mean(b_rev,'all') + dev*[-1 1];
                zlim(bounds), zticks(bounds), zlabel('Radians')
                zticklabels([sprintf("%.2fπ",bounds(1)) sprintf("%.2fπ",bounds(2))])
                title(sprintf("Phase detrending for shot %d",shot))

        end
    end

return