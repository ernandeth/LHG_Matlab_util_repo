clear, close all;
load spiralexampledata
% plot(kspacelocations);
[gdat kbkernal, kernaltable] = gridkb(kspacelocations,spiraldata,dcf,256,2,2);
[TF, ~, ~] = gridkb(kspacelocations,ones(size(spiraldata)),dcf,256,2,2);
im = fftshift(fft2(fftshift(gdat)));
im = abs(im)/500;
cmap = [0:127].'*[1 1 1] / 128;
colormap(cmap);
image(uint8(im));
colormap(cmap);
axis square;
