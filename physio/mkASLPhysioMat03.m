function PhysioMat = mkASLPhysioMat03(physfile, sampTime, disdaq, Nslices, TR, Ttag, taskM, AQ3D)
% function PhysioMat = mkASLPhysioMat(physfile, sampTime, disdaq, Nslices, TR, Ttag, taskM, [is3Dacquisition])
% This function generates a matrix of regressors to model the physiological
% effects on the ASL signal
%
% NOTE:  this version takes care of the slice timing by
% shifting the model for each slice.  It also decorrelates the 
% experimental design matrix from the physio matrix.
%
% physfile: names of the physio data file (4 column format)
%
% sampTime: sampling period for physio A/D
%
% disdaq: number of discarded images. these are not recorded in the Pfile, 
%    but the physio file contains data acquired during that time
%
% Nslices: number of slices in each frame
%
% TR: scanner image acquisition tims
%
% Ttag: time spent tagging (including the post inversion delay)
%
% taskM: task design matrix
%
% slTiming: flag to indicate timing on slices:  
%           0 = sequential, 
%           1 = 3D acquisition ... all are assumed to happen at center of slice
%           2 = assume that slice timing correction has been done and all slices are
%           acquired at the END of the TR.
%
% AQ3D: flag to indicate a 3D acquisition:  no slices!
%

if nargin < 8
    AQ3D=0;
end
AQ3D
phys = load(physfile);
Nimgs = round(size(phys,1)*sampTime / TR) - disdaq ;
Tslice = (TR-Ttag)/Nslices;
AQlength = (Nimgs+disdaq)*TR/sampTime

time = phys(1:AQlength,1);
resp = phys(1:AQlength,2);
card = phys(1:AQlength,3);

filtered_resp = smoothdata(resp, sampTime, 1, 9);
%plot(filtered_resp)

fprintf('\nExtracting respiratory phases ... finding resp peaks');
respPeaks = getPeaks(filtered_resp);

fprintf('\nMaking the resp phases vector');
rphases = zeros(size(resp));
rphases(respPeaks) = 1;
respPeakDist = diff(respPeaks);

for r=1:length(respPeaks)-1;
    rphases(respPeaks(r):respPeaks(r+1)) = [0:2*pi/respPeakDist(r): 2*pi];
end

fprintf('\nMaking the cardiac phases vector');
cardPeaks = find(card);
cphases = card;
cardPeakDist = diff(cardPeaks);

for r=1:length(cardPeaks)-1;
    cphases(cardPeaks(r):cardPeaks(r+1)) = [0:2*pi/cardPeakDist(r): 2*pi];
end

fprintf('\nMaking polynomial regressons (3rd order) for detrending');
t = (0:Nimgs-1)';
polynomial = [ones(size(t)) t/sum(t) t.^2/sum(t.^2)  t.^3/sum(t.^3) t.^4/sum(t.^4)];
polynomial = [ones(size(t)) t/sum(t) t.^2/sum(t.^2)  ];
mns = mean(polynomial,1); mns(1) = 0;
polynomial = polynomial - repmat(mns, length(t),1);

fprintf('\nRe-Sampling the phases at each slice')

%PhysioMat = zeros(Nslices, Nimgs, 14);
PhysioMat = zeros(Nslices, Nimgs, 12);


for slnum=0:Nslices -1
    switch slTiming
        case 0
            sl = slnum;        
        
        case 1
            sl = 1;  % this is when the center of Kspace is acquired

        case 2
            sl=Nslices;  % data are interpolated to the end of the TR
    end

    fprintf('\nresampling for slice ...%d',sl)
    ts = ((disdaq)*TR + Ttag + Tslice*sl)/sampTime : ...
        round(TR/sampTime): ...
        AQlength-1;
    ts = round(ts);

    cph = [sin(cphases(ts)) cos(cphases(ts)) sin(2*cphases(ts)) cos(2*cphases(ts)) ];
    rph = [sin(rphases(ts)) cos(rphases(ts)) sin(2*rphases(ts)) cos(2*rphases(ts)) ];
%     cph = [sin(cphases(ts)) cos(cphases(ts))  ];
%     rph = [sin(rphases(ts)) cos(rphases(ts))  ];
%     
    aslmod = ones(length(cph),1);
    aslmod(2:2:end, :) = -1;
%     DM = [polynomial cph rph  ];
%    DM = [polynomial cph rph  cph.*aslmod rph.*aslmod ];
    DM = [ polynomial cph rph   ];
    
    % decorrelate the matrix from the task

    DM = DM - taskM*pinv(taskM)*DM;
    
    DM = [aslmod DM];
    
    PhysioMat(slnum+1,:,:) = DM;
    
    
end
fprintf('\nCorrelation matrix');
rhos = corrcoef(DM)
subplot(211)
imagesc(squeeze(DM))
subplot(212)
imagesc(rhos); colorbar; title('correlation matrix for ASL physio matrix')
return
