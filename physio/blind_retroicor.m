function blind_retroicor(datafile, threshold, isASL)
% function blind_retroicor(datafile, threshold, isASL)
%
% this function identifies the brightest pixels (higher than threshold *SD)
% in a time series of images (datafile)as CSF or arteries, 
% and removes the average of those time courses from
% all the pixels in the image
%

% isASL = 1
%threshold =2;
discardRegs = 2;
P = getASLPhysMat(datafile,threshold);
if ~isASL
	P = [P(:,1) P(:,3)];
    discardRegs=1
end

%P(:,2) = P(:,2) - mean(P(:,2));
[m h] = read_img('mask');
mask = lightbox('mask');
figure; plot(P)
rmReg(datafile, P,discardRegs)
save PhysMat
return

