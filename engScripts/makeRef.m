function ref = makeRef(TR, duration, cond_duration, offset)
% function makeRef(TR, duration, condtion_duration, offset)
%
% this makes a very simple reference waveform for running correlation
% purposes and saves it as a file:  ReferenceWave.mat
%
%  warning This is  hard coded for ASL signals not BOLD !!
%

Nframes = duration/TR;


numTrials = duration/(2*cond_duration);
onsets = [offset: 2*cond_duration:duration-1]

ref = zeros(duration , 1);
for c=1:length(onsets)
    ref(onsets(c) +1 : (onsets(c) + cond_duration ) ) = 1;
end


% BOLD:
% h = spm_hrf(1, [6, 16, 1,1, 6,0,32]); 

% ASL:  shorter time to peak, little undershoot.
h = spm_hrf(1, [5, 16, 1,1, 20,0,32]); 
plot(h)
reg = conv(ref,h);
reg = reg(1:duration);


Xhires =  reg/max(reg);

% making the ASL matrix:
X = Xhires(1:TR:end, :);
X(1:2:end,:) = -X(1:2:end,:);

ref = -differencer(X, 4);
close

plot(X); hold on; plot(ref,'r')

save ReferenceWave.mat ref

return
