% load norm_job_100416.mat
% load norm_test_smooth.mat

rootDir = '/Users/mbu/Documents/fMRI_data/ASL_TMS_project/';
rootDir = '/home/data/asl/';
cd(rootDir)
load mfiles/norm_reslice100923.mat
 %{
        '090804jc'; % mean_sub - not flipped
        '090818jc'; % mean_sub - not flipped
    
        '090909jo'; % mean_sub is flipped; flipped = run flipper = 1
        '090917jo'; % mean_sub is flipped - suspect
        '091102jw'; % % mean_sub is flipped  - come back here too - missing rmean_sub
        '091111jw';% mean_sub is flipped

%}

subjs = [ ...
 
'091026bb';% mean_sub is flipped

'091104bb';% mean_sub is flipped
'091109ac'; % not flipped
'091118ac';% mean_sub is flipped
'091201el'; % mean_sub is flipped
'091207dh'; %  run 1 - not flipped, run 2 flipped? do it flipped - suspect
'091208jk';% not flipped  - suspect
'091130cl';  % not flipped
'091209cl';% not flipped
'091210el';% not flipped
'091216dh'; % mean_sub is flipped
'091217jk'; % mean_sub is flipped
'100125mh';  % - ok
'100203mh';  %  -ok
'100201ra';% - ok
'100208ep';% - ok

'100210ra';% - ok
 '100217ep';% - ok
'100308am';% - ok - missing a bit up top
'100317am';% - ok
'100322nt';% - ok
'100331nt';% - ok
'100426ab';% - ok - missing a bit up top
'100427as';% - suspect
'100505ab';
'100506as';
'100518as';
'100525kt';
'100525pc';
'100527as';
'100603kt';
'100603pc';
'100609hs';
'100615bk';
'100616nf';
'100618hs';
'100624bk';
'100625nf';
'100706lh';
'100720al';
'100726kn';
'100727vy';
'100729al';
'100804kn';
'100805vy';
'100817sh';
'100826sh';
];




doFlip=0;
doSmooth=0;

for s=1:size(subjs,1)
    
    subjDir = subjs(s,:)
    cd ([rootDir subjDir '/run_01']);
    %
    if doSmooth
        % smooth the images to help the normalization:
        in = ([rootDir  subjDir '/run_01/mean_sub.img']);
        out = ([rootDir  subjDir '/run_01/s_mean_sub.img']);
        spm_smooth(in,out,[ 8 8 8], 4 );
        
        in = ([rootDir  subjDir '/run_02/mean_sub.img']);
        out = ([rootDir  subjDir '/run_02/s_mean_sub.img']);
        spm_smooth(in,out,[ 8 8 8], 4 );
    end
    
    
    % make sure the headers are the same in all the files.
    
    
    
    h = read_hdr('mean_sub.hdr');
    write_hdr('mean_sub.hdr',h);
    h.datatype = 16;
    h.bits=32;
    h.dims=3;
    write_hdr('Zmap_0001.hdr',h);
    write_hdr('Zmap_0002.hdr',h);
    write_hdr('Zmap_0003.hdr',h);
    write_hdr('Zmap_0004.hdr',h);
    orig = h.origin;
    
    h = read_hdr('ConBhats.hdr');
    h.dims=4;
    h.tdim=4;
    h.origin = orig;
    write_hdr('ConBhats.hdr',h);
    write_hdr('ExpFlows.hdr',h);
    
    % flip images if necessary
    if doFlip
        zflipper('mean_sub' );
        zflipper('Zmap_0001');
        zflipper('Zmap_0002');
        zflipper('Zmap_0003');
        zflipper('Zmap_0004');
        
        zflipper('ConBhats')
        zflipper('ExpFlows');
    end

    % set up the normalization job and execute:
    
    jobs{1}.spatial{1}.normalise{1}.write.subj.matname{1} =    ....
        [ rootDir subjDir '/run_01/mean_sub_sn.mat'] ;
    
    jobs{1}.spatial{1}.normalise{1}.write.subj.resample = {
        [ rootDir subjDir '/run_01/mean_sub.img,1'] ,
        [ rootDir subjDir '/run_01/Zmap_0001.img,1'],
        [ rootDir subjDir '/run_01/Zmap_0002.img,1'],
        [ rootDir subjDir '/run_01/Zmap_0003.img,1'],
        [ rootDir subjDir '/run_01/Zmap_0004.img,1'],
        
        [ rootDir subjDir '/run_01/ConBhats.img,1'],
        [ rootDir subjDir '/run_01/ConBhats.img,2'],
        [ rootDir subjDir '/run_01/ConBhats.img,3'],
        [ rootDir subjDir '/run_01/ConBhats.img,4'],
        
        [ rootDir subjDir '/run_01/ExpFlows.img,1'],
        [ rootDir subjDir '/run_01/ExpFlows.img,2'],
        [ rootDir subjDir '/run_01/ExpFlows.img,3'],
        [ rootDir subjDir '/run_01/ExpFlows.img,4'],
        };
    
    spm_jobman('run', jobs);
    %
    
    if doFlip
        % make some links so that the names agree between those that needed
        % flippng and those that didn't
        ! ln -s wrExpFlows.img wExpFlows.img
        ! ln -s wrExpFlows.hdr wExpFlows.hdr
    end
    
    figure(55)
    [d h] = lightbox( [ rootDir subjDir '/run_01/mean_sub' ] );
    ov([],reshape(5* d, h.xdim, h.ydim, h.zdim),  32,32,6,0  );
    
    figure (56)
    
    [d h] = lightbox( [ rootDir subjDir '/run_01/wmean_sub'] );
    ov([],reshape(5* d, h.xdim, h.ydim, h.zdim),  26,32,25,0  );
    
    % make some links so that the names agree:
    
    
    %% run #2
    cd ([rootDir subjDir '/run_02']);
    %
    
    
    % make sure the headers are the same in all the files.
    
    
    
    h = read_hdr('mean_sub.hdr');
    write_hdr('mean_sub.hdr',h);
    h.datatype = 16;
    h.bits=32;
    h.dims=3;
    write_hdr('Zmap_0001.hdr',h);
    write_hdr('Zmap_0002.hdr',h);
    write_hdr('Zmap_0003.hdr',h);
    write_hdr('Zmap_0004.hdr',h);
    orig = h.origin;
    
    h = read_hdr('ConBhats.hdr');
    h.dims=4;
    h.tdim=4;
    h.origin = orig;
    write_hdr('ConBhats.hdr',h);
    write_hdr('ExpFlows.hdr',h);
    
    % flip images if necessary
    if doFlip
        zflipper('mean_sub' );
        zflipper('Zmap_0001');
        zflipper('Zmap_0002');
        zflipper('Zmap_0003');
        zflipper('Zmap_0004');
        
        zflipper('ConBhats')
        zflipper('ExpFlows');
    end
    
    
                
    jobs{1}.spatial{1}.normalise{1}.write.subj.matname{1} =    ....
        [ rootDir subjDir '/run_02/mean_sub_sn.mat'] ;
    
    
    jobs{1}.spatial{1}.normalise{1}.write.subj.resample = {
        [ rootDir subjDir '/run_02/mean_sub.img,1'] ,
        [ rootDir subjDir '/run_02/Zmap_0001.img,1'],
        [ rootDir subjDir '/run_02/Zmap_0002.img,1'],
        [ rootDir subjDir '/run_02/Zmap_0003.img,1'],
        [ rootDir subjDir '/run_02/Zmap_0004.img,1'],
        
        [ rootDir subjDir '/run_02/ConBhats.img,1'],
        [ rootDir subjDir '/run_02/ConBhats.img,2'],
        [ rootDir subjDir '/run_02/ConBhats.img,3'],
        [ rootDir subjDir '/run_02/ConBhats.img,4'],
        
        [ rootDir subjDir '/run_02/ExpFlows.img,1'],
        [ rootDir subjDir '/run_02/ExpFlows.img,2'],
        [ rootDir subjDir '/run_02/ExpFlows.img,3'],
        [ rootDir subjDir '/run_02/ExpFlows.img,4'],
        };
    
    spm_jobman('run', jobs);
    %
    
    % make some links so that the names agree:
    if doFlip
        ! ln -s wrExpFlows.img wExpFlows.img
        ! ln -s wrExpFlows.hdr wExpFlows.hdr
    end
    
    figure(55)
    [d h] = lightbox( [ rootDir subjDir '/run_02/mean_sub' ] );
    ov([],reshape(5* d, h.xdim, h.ydim, h.zdim),  32,32,6,0  );
    
    figure (56)
    
    [d h] = lightbox( [ rootDir subjDir '/run_02/wmean_sub'] );
    ov([],reshape(5* d, h.xdim, h.ydim, h.zdim),  26,32,25,0  );
    
    
    cd (rootDir);
end
