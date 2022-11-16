function aslsub_rel(imgname, M0frames) 

doDespike = 0;
doMask = 0;
msk = 1;

[mag h] = read_img(imgname);
        
    nonsel =mag(M0frames+1:2:end,:);
    sel = mag(M0frames+2:2:end,:) ;
    base = mag(1,:);
    
    if doDespike
        DespikeLevel = 2;
        c = sel;
        t = nonsel;
        
        mc = mean(c,2);
        stdmc = std(mc);
        mmc = mean(mc);
        badinds = find( abs(mc - mmc) > DespikeLevel *stdmc)
        c(badinds,:) = [];
        
        
        mt = mean(t,2);
        stdmt = std(mt);
        mmt = mean(mt);
        badinds = find( abs(mt - mmt) > DespikeLevel*stdmt)
        t(badinds,:) = [];
        
        sel= c;
        nonsel = t;
    end
    
    XRES = h.xdim;
    nonsel = reshape(mean(nonsel,1), XRES,XRES, h.zdim);
    sel = reshape(mean(sel,1), XRES,XRES, h.zdim);
    base = reshape(mean(base,1), XRES,XRES, h.zdim);

%     sel = smooth3(sel,'gaussian',[3 3 3]);
%     nonsel = smooth3(nonsel,'gaussian',[3 3 3]);
%     base = smooth3(base,'gaussian',[5 5 5]);
%     
    
    if doMask
        
            fprintf('making a mask')
            msk = base;
            msk(msk <= 1*mean(msk(:))) = 0;
            msk(msk>0) = 1;
            
            subplot(211)
            lightbox(msk); title('mask');
   
    end
    sel = sel ;
    nonsel = nonsel ;
    base = base ;
    
    d1 = 100*msk.*(nonsel-sel) ./ (base) ;
   
    % optional:  if you are not sure of the control-label order of the
    % images, force the difference to be positive on average
     if mean(d1(~isnan(d1))) < 0 
         d1 = -d1; 
     end

    tmph = h;
    tmph.tdim=1;
    write_img('mean_rsub.img', d1(:)*1000, tmph);

    subplot(212)
    lightbox(d1); axis ij;
    title('100 x [S_{nonsel} - S_{sel}] / (S_0)'); caxis([-2 2]);
    return
    
