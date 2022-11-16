warning off
rootDir = '/Users/mbu/Documents/fMRI_data/ASL_TMS_project/';
rootDir = '/home/data/asl/';
% % exclude:
%    234 data points?   X mat is 284...
%     '090804jc'; % mean_sub - not flipped
%     '090818jc'; % mean_sub - not flipped
%     '090909jo'; % mean_sub is flipped; flipped = run flipper = 1
%     '090917jo'; % mean_sub is flipped - suspect
    %'091102jw'; % % mean_sub is flipped  - come back here too - missing rmean_sub -> exclude
    %'091111jw';% mean_sub is flipped -> exclude
    
subjects = [ ...

    %'090804jc'; % pilot subject
    %'090909jo'; % different design matrix
    '091026bb';
    %'091102jw'; % over-performer -> did 5-back
    '091109ac';
    '091130cl'; % bad training data (exclusion debateable)
    '091201el';
    '091207dh';
    '091208jk';
    '100125mh'; % bad training data (exclusion debateable)
    '100201ra';
    '100208ep';
    '100308am';
    '100322nt'; % bad training data (exclusion debateable)
    %'100405ss' % scanner failure
    '100426ab';
    '100427as';
    '100518as';
    '100525kt';
    %'100518hb' % claustrophobic -> partial pre-test data and no post-data 
    '100525pc';
    %'100608jh' % tagging only on one side of the brain
    '100609hs';
    '100615bk';
    '100616nf';
    %'100706lh' % no show on post-test
    '100720al';
    '100726kn';
    '100727vy';
    '100817sh';
    '101109sw';
    
% post_test = [ ...
    %'090818jc'; % pilot subject
    %'090917jo'; % different design matrix
    '091104bb';
    %'091111jw'; % over-performer -> did 5-back
    '091118ac';
    '091209cl'; % bad training data (exclusion debateable)
    '091210el';
    '091216dh';
    '091217jk';
    '100203mh'; % bad training data (exclusion debateable)
    '100210ra';
    '100217ep';
    '100317am';
    '100331nt'; % bad training data (exclusion debateable)
    '100505ab';
    '100506as';
    '100527as';
    '100603kt';
    '100603pc';
    %'100617jh' % tagging only on one side of the brain
    '100618hs';
    '100624bk';
    '100625nf';
    '100729al';
    '100804kn';
    '100805vy';
    '100826sh';
    '101118sw';
    ];



Nsubs = size(subjects,1);
doNuisance=1;

alpha = load('/home/data/asl/mfiles/alphas2b.txt');
alpha = mean(alpha,2);

for mysub=1: Nsubs

    myalpha = alpha(mysub);
    
    if useSingleAlpha
        myalpha = 0.8;
    end
    
    for myrun=1:2
        close all
        fmri_name = subjects(mysub,:);
        
        [rootDir fmri_name filesep 'run_0' num2str(myrun)]
        cd([rootDir fmri_name filesep 'run_0' num2str(myrun)])
        
        
        %{
	if doNuisance
            close all
            [v h] = read_img('sub.img');
            vmap = var(v,0,1);
            h.tdim = 1;
            t = 2*median(vmap(:));
            vmap(vmap<t) = 0;
            write_img('variance.img', vmap,h);
            
            outname='wm_junk_PCA';
            
            ortho2005([],...
                'anat_file',[],...
                'tseries_file', 'sub.img', ...
                'ROIsize',1,...
                'spm_file', 'variance.img',...
                'spm_file2', [],...
                'output_name', outname, ...
                'doMovie',0 ...
                );
        end
        %}

        single_subj_fasl_script
        
        
    end
end
