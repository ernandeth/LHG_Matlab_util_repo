conditions = [1 2 3 4];
stimList = myDesign;
ISI = 1.2;
HRF = spm_hrf(ISI);

for i = 1:size(conditions,2)
   regressor = (stimList == conditions(i));
   regressor = conv(regressor, HRF);
   predMatrix(:,i) = regressor;
end
predMatrix = predMatrix(1:size(stimList),:);      	% eliminate extra values

