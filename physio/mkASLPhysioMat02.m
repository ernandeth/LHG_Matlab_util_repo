function PhysioMat = mkASLPhysioMat(physfile, sampTime, disdaq, Nslices, TR, Ttag, AQ3D)
% function PhysioMat02 = mkASLPhysioMat(physfile, sampTime, disdaq, Nslices, TR, Ttag, [is3Dacquisition])
% physfile: names of the physio data file (4 column format)
% sampTime: sampling period for physio A/D
% disdaq: number of discarded images. these are not recorded in the Pfile, but the physio file contains data
% acquired during that time
% Nslices: number of slices in each frame
% TR: scanner image acquisition tims
% Ttag: time spent tagging
% AQ3D: flag to indicate a 3D acquisition:  no slices!
%
% this version applies a simple differencing matrix to the physio data matrix 
%

if nargin < 7
    AQ3D=0;
end
AQ3D
phys = load(physfile);
Nimgs = round(size(phys,1)*sampTime / TR) - disdaq ;
Tslice = (TR-Ttag)/Nslices;
AQlength = (Nimgs+disdaq)*TR/sampTime

Dsimp = zeros(Nimgs/2 , Nimgs);
for c=1:Nimgs/2
    Dsimp(c, 2*c-1) = 1;
    Dsimp(c, 2*c) = -1;
end

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
PhysioMat = zeros(Nslices, Nimgs/2, 11);


for slnum=0:Nslices -1
    % 3D acquisitions don't have the slice timing issue
    % hard code this:
    if AQ3D
        sl = 1;  % this is when the center of Kspace is acquired
    else
        sl = slnum;
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
%    DM = [aslmod polynomial cph rph   ];
    
    DM = Dsimp * [cph rph ];
    DM = [polynomial(1:end/2,:)  DM];
    PhysioMat(slnum+1,:,:) =  DM;
    
    
end
fprintf('\nCorrelation matrix');
rhos = corrcoef(DM)
subplot(211)
imagesc(squeeze(DM))
subplot(212)
imagesc(rhos); colorbar; title('correlation matrix for ASL physio matrix')
return
