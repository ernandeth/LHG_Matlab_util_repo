function raw_corr = despike_vsasl3dflex(raw,M0frames,doRecycle,doSave)
    
    % Set defaults
    if nargin<1
        pfiles = dir('P*.7');
        raw = read_raw_3d(pfiles(1).name,0);
    end
    
    if nargin<2
        M0frames = 2;
    end
    
    if nargin<3
        doRecycle = 1;
    end
    
    if nargin<4
        doSave = doRecycle;
    end
    
    % If recycling, try to load despike data from file
    lastraw_corr = []; lastM0frames = 0;
    if doRecycle && isfile('./lastdespiked.mat')
        load ./lastdespiked.mat lastraw_corr lastM0frames
    end
    
    % Check that sizes and number of M0 frames match
    if (lastM0frames == M0frames) && (isequal(size(raw),size(lastraw_corr)))
        fprintf('\nSuccessfully loaded previously despiked data');
        raw_corr = lastraw_corr;
        doSave = 0;
    else
        fprintf('\nRejecting frame-wise outliers in k-space data');
        
        % Initialize corrected data array
        raw_corr = raw;
        
        % extract sizes/create frame arrays for control and label seperately
        nframes = size(raw,1);
        fv_con = (M0frames+1):2:nframes;
        fv_lab = (M0frames+2):2:nframes;
        ncoils = size(raw,3);
        
        % set up progress display
        if usejava('desktop')
            progbar = sprintf('\b|\n');
            fprintf([', coil-wise progress:\n' repmat('.',1,ncoils) ' (%d coils total)\n\t'], ncoils);
        else
            progbar = '';
            fprintf(' (cannot display coil-wise progress in nodesktop mode) ');
        end
        
        % initialize totalpts
        totalpts = 0;
        
        for coiln = 1:ncoils
            
            raw_corr_coil = squeeze(raw_corr(:,:,coiln));
            
            % find outliers for current coil
            [ir_con,ic_con] = find(isoutlier(abs(raw(fv_con,:,coiln)),1));
            [ir_lab,ic_lab] = find(isoutlier(abs(raw(fv_lab,:,coiln)),1));
            icu_con = unique(ic_con);
            icu_lab = unique(ic_lab);
            
            % interpolate over outliers for control frames
            for i = 1:length(icu_con)
                fvc = fv_con;
                col = icu_con(i);
                row = ir_con(ic_con==col);
                fvc(row) = [];
                
                raw_corr_coil(fv_con,col) = ...
                    interp1(fvc,raw(fvc,col,coiln),fv_con,'linear','extrap');
            end
            
            % ... and label frames
            for i = 1:length(icu_lab)
                fvc = fv_lab;
                col = icu_lab(i);
                row = ir_lab(ic_lab==col);
                fvc(row) = [];
                
                raw_corr_coil(fv_lab,col) = ...
                    interp1(fvc,raw(fvc,col,coiln),fv_lab,'linear','extrap');
            end
            
            totalpts = totalpts + length(icu_lab) + length(icu_con);
            raw_corr(:,:,coiln) = raw_corr_coil;
            fprintf('%s',progbar);
            
        end
        
        fprintf('\b (Done!)\n\t~%.2f%% of points removed', 100*totalpts/numel(abs(raw(:))));
    end
    
    if doSave
        lastraw_corr = raw_corr;
        lastM0frames = M0frames;
        save ./lastdespiked.mat lastraw_corr lastM0frames -v7.3
    end

end