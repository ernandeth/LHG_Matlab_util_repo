
rootDir = '/home/data/asl'

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
    '101109sw';
    '101118sw';
    ];

%handle only new subjects
subjects = [ ...
    '101207jb';
    '101213mz';
    '101216jb';
];

Nsubs = size(subjects,1);
cd(rootDir)
for mysub=1: Nsubs

    cd (subjects(mysub,:));
    cd('run_01')
    name = dir('vol_e*.nii');
    v = read_img(name(1).name);
    v = lightbox(reshape(v(1,:), 64,64,12));

    %ls *.dat
    
    str = input('is the top slice first?','s')
    if str=='y'
        sl = [12:-1:1];
        save sliceorder.dat sl -ascii
        !cp sliceorder.dat ../run_02
    end
    cd ../..
end


