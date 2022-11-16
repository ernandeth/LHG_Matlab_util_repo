warning off
rootDir = '/Users/mbu/Documents/fMRI_data/ASL_TMS_project/';
rootDir = '/home/data/asl/';


subjects = [ ...
    
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
'100608jh';
'100609hs';
'100615bk';
'100616nf'; %done
'100617jh';
'100618hs';
'100624bk';
'100625nf';
'100706lh';
'100720al'; %done
'100726kn'; %done
'100727vy'; %done
'100729al';
'100804kn';
'100805vy';
'100817sh'; %done
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

Nsubs = [];

for mysub=1:Nsubs
    for myrun=1:2
        
        fmri_name = subjects(mysub,:);
        
        [rootDir fmri_name filesep 'run_0' num2str(myrun)]
        cd([rootDir fmri_name filesep 'run_0' num2str(myrun)])        
        
%         !rm *.img.*
%         !rm  *Zmap* *Con* *Exp*
%         
        names = dir('*.img');
        for n=1:length(names)
            workFile =names(n).name
            [pth thename ext] = fileparts(workFile);
            zflipper(thename, 1);

        end
        sl = [12:-1:1];
        save sliceorder.dat sl -ascii
        !cp sliceorder.dat ../run_02
        
        ortho2005([],...
            'anat_file',[],...
            'tseries_file', 'sub.img', ...
            'ROIsize',1,...
            'spm_file', [],...
            'spm_file2', [],...
            'output_name', [], ...
            'doMovie',0 ...
            );
        str = pwd
        lightbox('mean_sub');
        title(str)
        pause
        close all
    end
end