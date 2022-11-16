% load norm_job_100416.mat
% load norm_test_smooth.mat

rootDir = '/home/data/asl/';
cd(rootDir)

pre_test = [ ...
    % '090804jc'; % mean_sub - not flipped
    % '090909jo'; % mean_sub is flipped; flipped = run flipper = 1
     '091026bb';% mean_sub is flipped
    % '091102jw'; % % mean_sub is flipped  - come back here too - missing rmean_sub
    '091109ac'; % not flipped
    %'091130cl';  % not flipped
    '091201el'; % mean_sub is flipped
    '091207dh'; %  run 1 - not flipped, run 2 flipped? do it flipped
    '091208jk';% not flipped
    '100125mh'; %- bad norm
    '100201ra';
    '100208ep';
    '100308am';
    '100322nt';
    '100426ab';
    '100427as';
    '100518as';
    '100525kt';
    '100525pc';
    '100609hs';
    '100615bk';
    '100616nf';
    '100720al';
    '100726kn';
    '100727vy';
    '100817sh';
    '101109sw';
    ];
post_test = [ ...
    % '090818jc'; % mean_sub - not flipped
    % '090917jo'; % mean_sub is flipped
     '091104bb';% mean_sub is flipped
    % '091111jw';% mean_sub is flipped
    '091118ac';% mean_sub is flipped
    %'091209cl';% not flipped
    '091210el';% not flipped
    '091216dh'; % mean_sub is flipped
    '091217jk'; % mean_sub is flipped
    '100203mh'; %- bad norm
    '100210ra';
    '100217ep';
    '100317am';
    '100331nt';
    '100505ab';
    '100506as';
    '100527as';
    '100603kt';
    '100603pc';
    '100618hs';
    '100624bk';
    '100625nf';
    '100729al';
    '100804kn';
     '100805vy';
    '100826sh';
    '101118sw';
    ];


underlay = '/home/data/asl/groupResults/rPET.nii';
for s=23:size(post_test,1)
    
    subjDir1=pre_test(s,:)
    subjDir2=post_test(s,:)
    
    Thr=2
    %{
 % check normlized maen subtractions from pretest and post test
    ortho2005([],...
        'ROItype', 'sphere',...
        'ROIsize', 15, ...
        'anat_file', [ rootDir subjDir1 '/run_01/mean_sub'], ...
        'spm_file', [ ], ...
        'threshold', Thr   ,...
        'spm_file2', [], ...
        'threshold2', Thr  ...
        );
    title(subjDir1);
    
    ortho2005([],...
        'ROItype', 'sphere',...
        'ROIsize', 15, ...
        'anat_file', [ rootDir subjDir2 '/run_01/mean_sub'], ...
        'spm_file', [ ], ...
        'threshold', Thr   ,...
        'spm_file2', [], ...
        'threshold2', Thr  ...
        );
    title(subjDir2);
    %}
    
    %  look at the stats maps for run 1
   
    
    ortho2005([],...
        'ROItype', 'sphere',...
        'ROIsize', 15,...
        'anat_file', underlay, ...
        'spm_file', [ rootDir subjDir1 '/run_01/wZmap_0004.img'], ...
        'threshold', Thr,...
        'spm_file2', [ rootDir subjDir2 '/run_01/wZmap_0004.img'], ...
        'threshold2', Thr);
    
    %% run #2
    

    ortho2005([],...
        'ROItype', 'sphere',...
        'ROIsize', 15, ...
        'anat_file', underlay, ...
        'spm_file', [ rootDir subjDir1 '/run_02/wZmap_0004.img'], ...
        'threshold', Thr   ,...
        'spm_file2', [ rootDir subjDir2 '/run_02/wZmap_0004.img'], ...
        'threshold2', Thr  ...
        );
    
    
    
    close all
    
end
