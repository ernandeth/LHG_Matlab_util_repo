function [phantom slice] = fakeBOLDslice( pSize, noiseLevel, tseries1, tseries2 )
% function [phantom slice]= fakeBOLDslice( pSize, noiseLevel, tseries1, tseries2 )
%
% create a slice of FMRI data with a square and a circular ROI
% the square gets the time course tseries1
% the circle gets the time course tseries2
% gaussian white noise is added to the whole thing
%

Nframes = length(tseries1);

slice= zeros(pSize);


% make ROI # 1

xcenter = round(pSize(1)/4);
ycenter = round(pSize(2)/4);
xwid = round(pSize(1)/8);
ywid = round(pSize(2)/8);

slice( xcenter-xwid:xcenter+xwid , ycenter-ywid:ycenter+ywid) = 1;

% make ROI #2

xcenter = round( 3*pSize(1)/4);
ycenter = round( 3*pSize(2)/4) ;
xwid = round(pSize(1)/8);
ywid = round(pSize(2)/8);

for x=xcenter-xwid:xcenter+xwid
    for y=ycenter-ywid:ycenter+ywid
    
        if sqrt((x-xcenter)^2 + (y-ycenter)^2 ) <= xwid,
            slice(x,y) = 2;
        end
    end
end

imagesc(slice);
slice = slice(:)';
Npix = length(slice);

phantom = repmat(slice, Nframes,1);

% stick in the time courses
for pix=1:Npix
    if slice(pix)==1
        phantom(:,pix) = tseries1;
    end
    if slice(pix)==2
        phantom(:,pix) = tseries2;
    end
end

% add noise
noise = sqrt(noiseLevel)*randn(Nframes, Npix);
phantom = phantom + noise;

return


