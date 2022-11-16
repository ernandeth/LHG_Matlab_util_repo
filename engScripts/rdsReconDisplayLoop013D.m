% This is called by rdsClientRecon during Real Time ASL FMRI
% the data, and the reconstruction params have already been loaded
% into the workspace.  This script gets called once per slice
%
% This file exists because it's easier to write out a separate
% file for matlab code than to make engEvalStrings
% within C programs using the Matlab Engine and recompile every time.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default ROI (i.e. when the program starts)
if InitialROIset==0
    x_roi=32;
    y_roi=32;
    z_roi=4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begin = (slicenum -1)* xres *xres +1;
finish = (slicenum)* xres *xres;

if DEBUG
    %fprintf('\n   Running in DEBUG mode... getting data from a variable named _raw_ ');
    %fprintf('\n   Make sure you set DEBUG=0 to operate in Real Time with the scanner');
    slice = raw(vnum,begin:finish);  % for testing!!
else
    % this is where the recon happens:
    % plot(abs(mxReconBuffer2)); title(['slice ' num2str(slicenum)]);drawnow
    slice = rt_sprec1(dat, args, scaninfo, kinfo);
    %fprintf('Reconned data! \n');
end

if SMOOTHFLAG
    % fast in-plane smoothing kernel using one neighbor on each side
    slice = reshape(slice, xres, xres);
    slice = slicesmooth(slice, 0.5);
    slice = reshape(slice, 1, xres*xres);
end

volumes(vnum,  begin:finish ) = slice(:)';
%fprintf('\n %d   %d ', vnum, slicenum);


if (slicenum==Nslices)

    if DO3D
        % this is where we do the FFT along the kz dimension if needed.
        tmp3d = reshape(volumes(vnum,:),xres*xres,Nslices);
        tmp3d = fftshift( fft( fftshift(tmp3d,2), [] , 2 ) );
        whos tmp3d volumes
        begin
        finish  volumes(vnum,  begin:finish ) = abs(tmp3d(:)); 
    end
    
    if DEBUG~=2
        subplot(221);
        lightbox(reshape(volumes(vnum,:),xres,xres,Nslices) );
        title(['raw data volume:' num2str(vnum)]);
    end
    %  save workspace;     %DEBUGGING ONLY!!!

    if (vnum==1)
        % masking by signal intensity and by the variance map
        slmask=volumes(1,:)>3*median(volumes(1,:));
        slmask = slmask .* varMask;
    end;

    if ~mod(vnum-3,6) & vnum>8

        cumasl = sum(allASL(vnum-5:vnum,:),1)/6;
        %cumasl = sum(allASL,1);

        if ~DOCORRELATION && ~DOREGRESSION
            ax2 = subplot(223);
            lightbox(reshape(cumasl, xres, xres, Nslices));
            title(['Cummulative ASL ...(six avgs)' num2str(vnum)]);

        end


    end




    if (vnum>2 && vnum<Nframes)
        % surround subtraction kernel: [1 -2 1]
        asl = volumes(vnum,:) - 2*volumes(vnum-1,:) + volumes(vnum-2,:) ;
        asl = asl*0.5*(-1)^vnum;
        asl = -asl.*slmask;

        allASL(vnum-2,:) = asl;


        mean_asl_course (vnum-2) = mean(asl(find(asl)));
        vol_temp=reshape(asl,xres,xres,Nslices);


        if DEMEANFLAG
            allASL(vnum-2,:) = asl - mean_asl_course (vnum-2) ;
        end

        if DEBUG~=2
            ax1 = subplot(222);
            lightbox(reshape(asl, xres, xres, Nslices),[-75 75],[]);
            title(['RT Surround Subtr ... (single avg)' num2str(vnum)]);
        end

        if DEBUG==2
            wmjunk=mean_asl_course;
        end
        
        if DOCORRELATION
            %%%
            running_corr
            %%%%
            
            if DEBUG~=2
                ax2 = subplot(223);
                act_lightbox( reshape(cumasl', xres,xres,Nslices), ...
                    reshape(corrMap',xres,xres,Nslices),...
                    [0 75], [0.4 1],[]);
                title(['Updated Correlation and Mean Subtraction ...' num2str(vnum)]);
            end
        end


        if DEBUG~=2
            % Find the coordinates of where I clicked
            myclick = get(gcf,'CurrentPoint');
            figsize = get(gcf,'Position');
            if myclick(2) > figsize(4)/2
                myAxis = ax1;
            else
                myAxis = ax2;
            end

            myclick = get(myAxis,'CurrentPoint');
            xc = abs(myclick(1,1));
            yc = abs(myclick(1,2));
            myrow = ceil(yc/xres) ;
            mycol = ceil(xc/xres) ;

            x_roi = round(rem(abs(xc),xres));
            y_roi = round(rem(abs(yc),xres));
            z_roi = Ncols * (myrow-1) + mycol;

            %        [myrow mycol x_roi y_roi z_roi]


            InitialROIset=1;

            if x_roi < 2, x_roi=2; end
            if y_roi < 2, y_roi=2; end
            if z_roi < 1, z_roi=1; end
            if x_roi > xres-2,  x_roi=xres-1; end
            if y_roi > xres-2,  y_roi=xres-1; end
            if z_roi > Nslices, z_roi=Nslices; end

            %
            [xx,yy,zz] = ndgrid([x_roi-1:x_roi+1], [y_roi-1 : y_roi+1], z_roi*[ 1 1 1]);
            nlist = length(xx(:));
            fx = reshape(xx,nlist,1);
            fy = reshape(yy,nlist,1);
            fz = reshape(zz,nlist,1);

            inds = sub2ind( [xres, xres, Nslices], fx, fy, fz);
            inds = inds(inds>0 & inds <= Npix);

            if myAxis == ax1
                wm_inds = inds;
            else
                ref_inds = inds;
            end

            wmjunk = mean(allASL(:,wm_inds),2);
            wmjunk = wmjunk - mean(wmjunk);

            mean_ROI_course = mean(clean(:,ref_inds),2);
            %%

            subplot(224)
            plot(mean_asl_course(1:vnum-2));
            hold on;
            plot(mean_ROI_course(1:vnum-2),'r');

            plot(wmjunk(1:vnum-2),'g');


            axis([1 Nframes -50 50]);
            hold off
            title('the mean ASL time course so far ...');

        end


    end

    
    if length(wmjunk)>1 && DOREGRESSION
        %%
        tic
        running_regression
        tcrunch(vnum) = toc;
        %%
        
        ax2 = subplot(223);

        act_lightbox( reshape(cumasl', xres,xres,Nslices), ...
            reshape(corrMap',xres,xres,Nslices),...
            [-25 25], [3 10],[]);
        title(['Updated T scores and Mean Subtraction ...' num2str(vnum)]);

        subplot(224)
        hold on
        plot(10*ref1,'k');
        hold off
    end


    if SAVEMOVIE
        movavi = addframe(movavi, getframe(gcf));
    end

    drawnow
    vnum = vnum+1;
    slicenum = 0;


end


slicenum = slicenum +1;

% make some screeenshots:
if vnum==(Nframes/2) && slicenum==1 && DEBUG
    print -depsc halfway
end
% make some screeenshots:
if vnum==(Nframes) && slicenum==1 && DEBUG
    print -depsc last
    save Tmap corrMap
end


if vnum==(Nframes) && slicenum==1

    if ~DEBUG
        title(sprintf('\n %d frames acquired.  saving RTdata.mat file \n', Nframes));
        save RealTimeASL.mat volumes allASL
    end

    if SAVEMOVIE
        movavi = close(movavi)
    end



end


