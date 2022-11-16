TR = 5;
hrf = spm_hrf(TR);
task = zeros(160,1);
onsets = ([25 60 200 240 400 440 600 640 780] - 20  )/TR;
durations = [35 140 40 160 40 160 40 140 40] / TR;

for c = 2:4:length(onsets)
	task(onsets(c):onsets(c)+durations(c))=1;
end

aslmod = ones(size(task));
aslmod(2:2:end) = -1;

distractor = task;

task = conv(task,hrf);
task = task(1:length(aslmod));
btask = task;

task = task .* aslmod;

distractor(:) =0;
for c = 4:4:length(onsets)
	distractor(onsets(c):onsets(c)+durations(c))=1;
end
distractor = conv(distractor,hrf);

bdistractor = distractor(1:length(aslmod)) ;
distractor = distractor(1:length(aslmod)) .* aslmod;

fixation= task;
fixation(:) =0;
for c = 1:2:length(onsets)
	fixation(onsets(c):onsets(c)+durations(c))=1;
end
fixation = conv(fixation,hrf);

bfixation = fixation(1:length(aslmod)) ;
fixation = fixation(1:length(aslmod)) .* aslmod;

baseline = aslmod;

save runCD_TR5_fixation.txt fixation -ascii
save runCD_TR5_bfixation.txt bfixation -ascii
save runCD_TR5_baseline.txt baseline -ascii
save runCD_TR5_task.txt task -ascii
save runCD_TR5_distractor.txt distractor -ascii
save runCD_TR5_btask.txt btask -ascii
save runCD_TR5_bdistractor.txt bdistractor -ascii

DM = [fixation bfixation baseline  distractor  task btask bdistractor];          
save DesMatCD.txt DM -ascii
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
TR = 5;
hrf = spm_hrf(TR);
task = zeros(160,1);
onsets = ([25 60 220 260 400 440 600 620 780] - 20  )/TR;
durations = [35 160 40 140 40 140 40 160 40] / TR;

for c = 2:4:length(onsets)
	task(onsets(c):onsets(c)+durations(c))=1;
end

aslmod = ones(size(task));
aslmod(2:2:end) = -1;

distractor = task;

task = conv(task,hrf);
task = task(1:length(aslmod));
btask = task;

task = task .* aslmod;

distractor(:) =0;
for c = 4:4:length(onsets)
	distractor(onsets(c):onsets(c)+durations(c))=1;
end
distractor = conv(distractor,hrf);

bdistractor = distractor(1:length(aslmod)) ;
distractor = distractor(1:length(aslmod)) .* aslmod;

fixation= task;
fixation(:) =0;
for c = 1:2:length(onsets)
	fixation(onsets(c):onsets(c)+durations(c))=1;
end
fixation = conv(fixation,hrf);

bfixation = fixation(1:length(aslmod)) ;
fixation = fixation(1:length(aslmod)) .* aslmod;

baseline = aslmod;

save runAB_TR5_fixation.txt fixation -ascii
save runAB_TR5_bfixation.txt bfixation -ascii
save runAB_TR5_baseline.txt baseline -ascii
save runAB_TR5_task.txt task -ascii
save runAB_TR5_distractor.txt distractor -ascii
save runAB_TR5_btask.txt btask -ascii
save runAB_TR5_bdistractor.txt bdistractor -ascii

DM = [ fixation bfixation baseline  distractor  task btask bdistractor];
save DesMatAB.txt DM -ascii

