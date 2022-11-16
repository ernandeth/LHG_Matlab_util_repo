% use this for realignment:
% ! mcflirt -in arun_05.nii  -out rarun_05.nii -refvol 100  -cost normcorr -verbose 1 -stats -plots -mats
%
% smoother('rarun_05.nii', 2);
%

TR = 1;
% notes - we collected 300 seonds of data but the paradigm ended after 255.
% clipped off the last 45 ffrom the end of the rarun_01.nii file

duration = 310 * TR;
fixtime= [...
	9 9 9 9 9 9 ...
	10 10 10 10 10 10  ...
	5 5 5 5 5 5 ...
	7 7 7 7 7 7 ...
	6 6 6 6 6 6 ...
	8 8 8 8 8 8 ];

acttime = ones(size(fixtime));

D = zeros(duration / TR , 2);
t = 0;
for c=1:length(fixtime)
    start = t + fixtime(c)/TR
    stop = t + fixtime(c)/TR + acttime(c)/TR
	t=stop
	D(start : stop, 2) = 1;
end

D(:,1) = 1;
D = D(1:duration,:);

h = spm_hrf(TR);
reg = D(:,2);
reg = conv(reg,h);
reg = reg(1:duration/TR);
reg = reg - mean(reg);

D(:,2) =  reg;

spmJr('srarun_05',D,[0 1]);

