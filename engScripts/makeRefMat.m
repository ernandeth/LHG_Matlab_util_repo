TR = 0.9;
TR=0.7;
TR=4;

onsets{1} = [20:50:320];
durations{1} = 20*ones(size(onsets{1}));

ref = buildDesMat(TR, 320, onsets, durations, 1);

ref = differencer(ref,4);

ref = -ref(:,4);
plot(ref);

save ReferenceWave.mat ref

