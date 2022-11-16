function grad = dotrap(area,gslew,gamp,gts)
%dotrap - generates a trap waveform
% units must match
% area = mT/m * ms
% gslew = mT/m/ms
% gamp = mT/m
% gts = ms

% determine if trap or tri
maxrampl = ceil(gamp/gslew/gts);
maxtriarea = maxrampl * gamp * gts;
if (area < maxtriarea) 
    rampl = ceil(sqrt(abs(area)/gslew)/gts);
    shape = [(0:rampl) ((rampl-1):-1:0)];
    grad = shape/sum(shape*gts)*area;
else
    flatl = ceil((abs(area)-maxtriarea)/gamp/gts);
    shape = [(0:maxrampl)/maxrampl ones([1 flatl]) ((maxrampl-1):-1:0)/maxrampl];
    grad = shape/sum(shape*gts)*area;
end
end

