TR=4;
lag=0;
exp_duration=112*4;

%% Building deisgn matrix before anything else
onsets{1} = [6: 56: 448] - lag  ;  % Craving
durations{1} = 20*ones(size(onsets{1}));

onsets{2} = [34: 56: 448]-lag ;
durations{2} = 20*ones(size(onsets{1}));


X = buildDesMat(TR, exp_duration, onsets, durations, 1);
DX = -differencer(X,4);

DX(:,1) = [];
imagesc(DX);
%%


Pfile='P09728.7';


args.inFile = Pfile;
args.doDespike = 0;
args.doRecon=2;  % do 3D recon
args.doSliceTime=0;

args.doRealign = 0;
args.smoothSize= 1;  % do mrfilter
args.subType = 2;
args.physCorr = 2;
args.physFile='';

args.doGLM = 1;
args.designMat = DX; % load('schizTMS_designmat.dat');
args.isSubtracted = 1;  % this means don't do any further subtraction
args.contrasts = [0 0 1 0 0;   0 0 0 -1 1];
args.doQuant = 1;
args.doGlobalMean = 0;

args.aslParms.TR = 4;
args.aslParms.Ttag = 1.8;
args.aslParms.Tdelay = 1.4;
args.aslParms.Ttransit = 1.2;
args.aslParms.inv_alpha = 0.8;
args.aslParms.disdaqs = 0;

args.doLightbox = 1;
args.doOrtho = 1;

args.doFlip = 0;
args.subOrder = 1;

% new stuff:
args.doBuildAnatomy =1;
args.anatomyDir = '../strucs/';
args.doCoreg = 1;
args.overlayfile='t1overlay.nii';
args.spgrfile='t1spgr.nii';
args.doTransfer=1;

asl_spm02 %(args)
