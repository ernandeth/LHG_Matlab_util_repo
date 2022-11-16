% this is a batch file to flip the subjects along teh z direction:

rootDir = '/Users/mbu/Documents/fMRI data/ASL TMS project/';

subjs = [ ...
%         %'090804jc'; % mean_sub - not flipped
%         %'090818jc'; % mean_sub - not flipped
%         '090909jo'; % mean_sub is flipped; flipped = run flipper = 1
%         '090917jo';  % mean_sub is flipped
%         '091026bb'; % mean_sub is flipped
%         '091102jw';  % mean_sub is flipped  - come back here too - missing rmean_sub
%         '091104bb';  % mean_sub is flipped
%         %'091109ac'; % not flipped
%         '091111jw'; % mean_sub is flipped
%         '091118ac'; % mean_sub is flipped
%         %'091130cl';  % not flipped
%         '091201el'; % mean_sub is flipped
%         %'091207dh'; %  run 1 - not flipped, run 2 flipped? do it flipped
%         '091208jk'; % not flipped
%         %'091209cl'; % not flipped
%         %'091210el'; % not flipped
%         '091216dh';  % mean_sub is flipped
%         '091217jk';  % mean_sub is flipped
        
        '091026bb';% mean_sub is flipped
%  '091102jw'; % % mean_sub is flipped  - come back here too - missing rmean_sub
% '091111jw';% mean_sub is flipped
'091104bb';% mean_sub is flipped

    ];


for s=1:size(subjs,1)
    
    subjDir = subjs(s,:)
    cd ([rootDir subjDir '/run_01/'])
    zflipper('ExpFlows');
    zflipper('ExpFlow_vars');
    !mv rExpFlows.img ExpFlows.img
    !mv rExpFlows.hdr ExpFlows.hdr
    !mv rExpFlow_vars.img ExpFlow_vars.img
    !mv rExpFlow_vars.hdr ExpFlow_vars.hdr

end
