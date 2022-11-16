%function pairedTtest_spmJr
% rootDir = '/net/martin/data/groupResults/nback/'
rootDir = '/Users/mbu/Documents/fMRI_data/ASL_TMS_project/';
rootDir = '/home/data/asl/';

doTest = 0;
doFigures = 1;
useConfound=1;

pre_test = [ ...
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
    '101207jb';
    '101213mz';
    ];
post_test = [ ...
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
    '101216jb';
    '101222mz';
    ];


destinationDir =[ rootDir filesep '4back_1back']
FLC = 1;
close all
for mycontrast = 2:3;

    switch (mycontrast)
        case 1
            destinationDir = [ rootDir filesep 'groupResults/baseline']
            FLC = 1;  % first level contrast:   1 = baseline,     4 = 4back-1back

        case 2
            destinationDir = [ rootDir filesep 'groupResults/4back']
            FLC = 2;  % first level contrast:   1 = baseline,     4 = 4back-1back

        case 3
            destinationDir = [ rootDir filesep 'groupResults/1back']
            FLC = 3;  % first level contrast:   1 = baseline,     4 = 4back-1back

        case 4
            destinationDir = [ rootDir filesep 'groupResults/4back_1back']
            FLC = 4;  % first level contrast:   1 = baseline,     4 = 4back-1back
    end

    if ~exist(destinationDir)
        eval(['!mkdir' destinationDir])
    end

    if doTest
        Nsubjects=size(pre_test,1);
        all_subjs = [];
        figure;
        % load all the flow estimates from the Nth (FLC) contrast:
        for s=1:Nsubjects

            %smoother3(  [ rootDir   pre_test(s,:) '/run_01/wExpFlows.img'], 5);

            flows = read_img(   [ rootDir   pre_test(s,:) '/run_01/swExpFlows']);
            all_subjs =  [all_subjs ; flows(FLC,:)];


            %smoother3(  [ rootDir   pre_test(s,:) '/run_02/wExpFlows.img'], 5);

            [flows hdr] = read_img(   [ rootDir   pre_test(s,:) '/run_02/swExpFlows']);
            all_subjs =  [all_subjs ;  flows(FLC,:)];

            %lightbox(reshape(flows(FLC,:),hdr.xdim, hdr.ydim, hdr.zdim)); drawnow
        end

        for s=1:Nsubjects

           %smoother3(   [ rootDir   post_test(s,:) '/run_01/wExpFlows.img'],5);

            [flows hdr] = read_img([ rootDir   post_test(s,:) '/run_01/swExpFlows']);
            all_subjs =  [all_subjs ; flows(FLC,:)];

           %smoother3( [ rootDir   post_test(s,:) '/run_02/wExpFlows.img'],5);

            [flows hdr] = read_img( [ rootDir   post_test(s,:) '/run_02/swExpFlows']);
            all_subjs =  [all_subjs ; flows(FLC,:)];

            % lightbox(reshape(flows(FLC,:),hdr.xdim, hdr.ydim, hdr.zdim)); drawnow
        end

        meanFlows = zeros(Nsubjects*4,1);
        for p=1:Nsubjects*4
            tmp = all_subjs(p,:);
            tmp = tmp(isfinite(tmp));
            tmp = tmp(find(tmp));
            meanFlows(p)= mean(tmp);
        end
        meanFlows = meanFlows-mean(meanFlows);
        meanFlows = meanFlows /max(meanFlows);

        % now build the design matrix
        X = [1  1]';
        X = kron(eye(Nsubjects),X);
        X = [X ; X];
        r = zeros(Nsubjects*4,2);
        r(1:end/2, 1) = 1;
        r(end/2+1:end, 2) = 1;
        X = [X r];

        % include inversion efficiency as
        % a regressor
        % read in the raw alphas from martin
        alphareg = load('/home/data/asl/mfiles/alphaRegressor.txt');
        %performance = load('/home/data/asl/mfiles/performRegressor.txt');

        C = [...
            zeros(1, Nsubjects) 1 -1   ;  % pre - post
            zeros(1, Nsubjects) -1 1   ;  % post - pre  (redundant)
            zeros(1, Nsubjects)  1 0   ;  % pre-test only
            zeros(1, Nsubjects)  0 1   ;  % post-test only
            ];

        if useConfound
            %X = [X alphareg performance];
            X = [X alphareg ];
            C = [...
                zeros(1, Nsubjects) 1 -1  0  ;  % pre - post
                zeros(1, Nsubjects) -1 1  0 ;  % post - pre  (redundant)
                zeros(1, Nsubjects)  1 0  0 ;  % pre-test only
                zeros(1, Nsubjects)  0 1  0 ;  % post-test only
                zeros(1, Nsubjects)  0 0 1  ;  % confound only
                % zeros(1, Nsubjects)  0 0 0  1;  % confound only
                ];
        end
        
        imagesc(X)

        % now do the analysis proper
        hdr.tdim = 1;
        flags.doWhiten = 0;
        flags.header = hdr;

        cd (destinationDir)
        % clean up first
        !rm *.img *.hdr

        % let's write the perfusion effects from all into a file so we can look at it later
        h = hdr;
        h.tdim = Nsubjects*4;
        write_img('groupFlowEffects.img', all_subjs, h);

        % sanity check: there should be no difference between run 1 and 2
        % in the baseline
        saneCon = zeros(1, size(X,2));
        saneCon(1:2:Nsubjects)=1;
        saneCon(2:2:Nsubjects)=-1;

        spmJr_flex(all_subjs, X, C, flags);


        lightbox('Zmap_0002',[-3 3],[]);
    end


    if doFigures

        cd(destinationDir)

        % make a pretty figure:
        underlay = '/home/data/asl/groupResults/rPET.nii';
        underlay = '/home/data/asl/groupResults/baseline/ConBhat_0003';
        

        %[f h] = read_img([rootDir '091216dh/run_01/wExpFlows']
        [f h] = read_img(underlay);
        f0 = reshape(f, h.xdim, h.ydim, h.zdim);


        feffect = read_img2([ destinationDir '/ConBhat_0004']);
        msk = read_img2([ destinationDir '/Zmap_0004.img']);

        figure('Position', [1 200*(mycontrast-1) 300 300])
        subplot(211); hist(msk(find(msk)),100);
        subplot(212); lightbox(msk, [-5 5], 6);

        msk(abs(msk)<2.5) = 0;
        msk(abs(msk)>=2.5) = 1;

        % Keep only the grey mater:
        msk2 = read_img2(underlay);
        msk2(abs(msk2)<25) = 0;
        msk2(abs(msk2)>=25) = 1;

        feffect = feffect.*msk.*msk2;
        %f0 = f0.*msk2;
        %feffect = feffect.*msk2;

        slices = [ 24: 2: 36]
        f0_b = f0(:,:,slices);
        feffect_b = feffect(:,:,slices);

        figure('Position', [100 200*(mycontrast-1) 1500 300])

        act_lightbox(f0_b, feffect_b, [0 80], [0.01  5], 1);
        title(destinationDir)
        print -dtiff post_flowsZcrit=2.5
    end

end

