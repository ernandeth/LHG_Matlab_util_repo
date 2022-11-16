% This is called by rdsClientRecon during Real Time ASL FMRI
% the data, and the reconstruction params have already been loaded
% into the workspace.  This script gets called once per slice
%
% This file exists because it's easier to write out a separate
% file for matlab code than to make engEvalStrings
% within C programs using the Matlab Engine and recompile every time.

warning off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default ROI (i.e. when the program starts)
if InitialROIset==0
    x_roi=32;
    y_roi=32;
    z_roi=4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



rbeg = (slicenum -1)* 2*scaninfo.ndat +1;
rfin = (slicenum)* 2*scaninfo.ndat;

if DEBUG
    %fprintf('\n   Running in DEBUG mode... getting data from a variable named _raw_ ');
    %fprintf('\n   Make sure you set DEBUG=0 to operate in Real Time with the scanner');
    dat = allRawData(vnum,rbeg:rfin);  % for testing!!
    
else
    % this is where the recon happens:
    % plot(abs(mxReconBuffer2)); title(['slice ' num2str(slicenum)]);drawnow
    allRawData(vnum, rbeg:rfin) = dat;
end





if (slicenum==Nslices)
    
    if DORECON
        % Try to do the spiral recon in parallel.
        parfor sl=1:Nslices
            slice = rt_sprec1(dat, args, scaninfo, kinfo);
            
            %fprintf('rt_sprec1 .... Reconned data! \n');
            
            if SMOOTHFLAG
                % fast in-plane smoothing kernel using one neighbor on each side
                % slice = reshape(slice, xres, xres);
                slice = slicesmooth(slice, 0.5);
                %slice = reshape(slice, 1, xres*xres);
            end
            
            if DO3D
                % the order of the kspace pancakes is not sequential.  we start at kz=0
                % and go out in k-space.
                tmp3d(:,:,order(slicenum)) = slice;
                % order(slicenum)
            else           
                slice = slice(:)';
                
                begin = (sl -1)* xres *xres +1;
                finish = (sl)* xres *xres;
                volumes(vnum,  begin:finish ) = abs(slice);
            end
            %fprintf('\n %d   %d ', vnum, slicenum);
        end
        
    end
    
    if DO3D && DORECON
        % this is where we do the FFT along the kz dimension if needed.
        tmp3drec = fftshift(abs(fft(fftshift(tmp3d,3) , [] , 3) ),3 );
        volumes(vnum, : ) = tmp3drec(:);
    end
    
    if DEBUG~=2
        % update the raw data display figure here
        figure(1);
        lightbox(reshape(volumes(vnum-1,:),xres,xres,Nslices),[],3 );
        title(['raw data volume:' num2str(vnum-1)]);
        drawnow
    end
    
    
    %  save workspace;     %DEBUGGING ONLY!!!
    
    if (vnum==1)
        % masking by signal intensity and by the variance map
        slmask=volumes(1,:)>3*median(volumes(1,:));
        slmask = slmask .* varMask;
    end;
    
    if ~mod(vnum-3,6) & vnum>8
        
        cumasl = sum(allASL,1)/(vnum-2);
        
        if ~DOCORRELATION && ~DOREGRESSION
            ax2 = figure(3);
            lightbox(reshape(cumasl, xres, xres, Nslices),[-200 200],[3]);
            title(['Cummulative ASL ...(six avgs)' num2str(vnum)]);
            drawnow
        end
    end
    
    if (vnum>2 && vnum<=Nframes)
        
        
        % surround subtraction kernel: [1 -2 1]
        asl = volumes(vnum,:) - 2*volumes(vnum-1,:) + volumes(vnum-2,:) ;
        asl = asl*0.5*(-1)^vnum;
        
        asl = asl.*slmask;
        if DEBUG~=2
            
            ax1 = figure(2);
            lightbox(reshape(asl, xres, xres, Nslices),[-200 200],3);
            title(['RT Surround Subtr ... (single avg)' num2str(vnum)]);
            drawnow
            
        end
        
        allASL(vnum-2,:) = asl;
        
        mean_asl_course (vnum-2) = mean(asl(find(asl)));
        
        if DEMEANFLAG
            
            allASL(vnum-2,:) = asl - mean_asl_course (vnum-2) ;
        end
        
        
        if DEBUG==2
            wmjunk=mean_asl_course;
        end
        
        if DEBUG~=2
            % Find the coordinates of where I clicked
            
            figure(2)
            myclick = get(gca,'CurrentPoint');
            xc = abs(myclick(1,1));
            yc = abs(myclick(1,2));
            hold on, plot(xc, yc,'g*'), hold off
            [wm_inds] = rdsfindcoords(xc, yc, xres, xres, Nslices);
            
            
            figure(3)
            myclick = get(gca,'CurrentPoint');
            xc = abs(myclick(1,1));
            yc = abs(myclick(1,2));
            hold on, plot(xc, yc,'g*'), hold off
            [ref_inds] = rdsfindcoords(xc, yc, xres, xres, Nslices);
            
            
            figure(4)
            % this is the nuisance time course from the asl map
            wmjunk = mean(allASL(:,wm_inds),2);
            junkreg = wmjunk(1:vnum-2);
            junkreg = junkreg - mean(junkreg);
            junkreg = junkreg / (eps+max(junkreg));
            plot(10*junkreg,'g');
            hold on;
            
            % plot the mean in the ROI of the cumulative map
            mean_ROI_course = mean(clean(:,ref_inds),2);
            mean_ROI_course = mean_ROI_course(1:vnum-2);
            mean_ROI_course = mean_ROI_course-mean(mean_ROI_course);
            mean_ROI_course = mean_ROI_course / (eps+max(mean_ROI_course));
            
            plot(10*mean_ROI_course,'r');
            
            plot(mean_asl_course(1:vnum-2),'b');
            
            axis([1 Nframes -50 50]);
            hold off
            
        end
        
        
    end
    
    
    if length(wmjunk)>1 && DOREGRESSION && vnum>5
        %%
        
        tic
        running_regression
        tcrunch(vnum) = toc;
        %%
        
        ax2 = figure(3);
        
        act_lightbox( reshape(cumasl', xres,xres,Nslices), ...
            reshape(corrMap',xres,xres,Nslices),...
            [-50 50], [2.5 20],3);
        title(['Updated T scores and Mean Subtraction ...' num2str(vnum)]);
        
        figure(4)
        hold on
        plot([nn:N], 10*ref1,'k');
        hold off
        drawnow
        
    end
    
    if DOCORRELATION
        %%%
        
        running_corr
        %%%%
        
        if DEBUG~=2
            
            ax2 = figure(3);
            act_lightbox( reshape(cumasl', xres,xres,Nslices), ...
                reshape(corrMap',xres,xres,Nslices),...
                [5 50], [0.4 1],3);
            
            title(['Updated Correlation and Mean Subtraction ...' num2str(vnum)]);
            
            
            
            figure(4)
            hold on
            plot(10*ref1,'k');
            hold off
            drawnow
        end
    end
    
    
    figure(3)
    if SAVEMOVIE
        movavi = addframe(movavi, getframe(3));
    end
    
    
    % make some screeenshots:
    if vnum==(Nframes/2)  && DEBUG
        
        print -depsc halfway
    end
    % make some screeenshots:
    if vnum==(Nframes) && slicenum==1 && DEBUG
        
        print -depsc last
        save Tmap corrMap
        matlabpool close
    end
    
    
    if vnum==(Nframes)
        
        if ~DEBUG
            title(sprintf('\n %d frames acquired.  saving RTdata.mat file \n', Nframes));
            save RealTimeASL.mat volumes allASL allRawData;
        end
        
        if SAVEMOVIE
            movavi = close(movavi)
        end
        
        
        
    end
    
    
    % reset the slice number index and increment the volume counter
    slicenum=1;
    vnum=vnum+1;
    
else
    slicenum = slicenum +1;
end


