% This is called by rdsClientRecon during Real Time ASL FMRI
% the data, and the reconstruction params have already been loaded
% into the workspace.  This script gets called once per slice
%
% This file exists because it's easier to write out a separate
% file for matlab code than to make engEvalStrings
% within C programs using the Matlab Engine and recompile every time.

DEBUG = 0;

begin = (slicenum -1)* xres *xres +1;
finish = (slicenum)* xres *xres;

if DEBUG
    %fprintf('\n   Running in DEBUG mode... getting data from a variable named _raw_ ');
    %fprintf('\n   Make sure you set DEBUG=0 to operate in Real Time with the scanner');
%    slice = raw(vnum,begin:finish);  % for testing!!
else
    % this is where the recon happens:
    slice = rt_sprec1(dat, args, scaninfo, kinfo);
end



volumes(vnum,  begin:finish ) = slice(:)';
%fprintf('\n %d   %d ', vnum, slicenum);

  
 if (slicenum==Nslices)

    subplot (221)
    lightbox(reshape(volumes(vnum,:),xres,xres,Nslices) );
    title(['raw data volume:' num2str(vnum)]);
    
 

%    save workspace;
if (vnum==1)
    slmask=volumes(1,:)>3*median(volumes(1,:));
end;




    if (vnum>2 && vnum<Nframes)
        asl = volumes(vnum,:) - 2*volumes(vnum-1,:) + volumes(vnum-2,:) ;
        asl = asl*0.5*(-1)^vnum;
        asl = asl.*slmask;
        
   allASL(vnum-2,:) = asl;
       mean_asl_course (vnum-2) = mean(asl);
       
        subplot(222)
        lightbox(reshape(asl, xres, xres, Nslices),[-50 100],4);
        title(['RT Surround Subtr ...' num2str(vnum)]);
        
        if DEBUG, pause(0.1), end
            
    end

    if ~mod(vnum-3,6) & vnum>8

        cumasl = sum(allASL(vnum-5:vnum,:),1);
        %cumasl = sum(allASL,1);
        
        subplot(223)
        lightbox(reshape(cumasl, xres, xres, Nslices));
        title(['Cummlative ASL ...' num2str(vnum)]);



        subplot(224)
        plot(mean_asl_course(1:vnum-2));
        title('the mean ASL time course so far ...');
        drawnow

        
    end
        vnum = vnum+1;
     slicenum = 0;
    
 end


slicenum = slicenum +1;

