function hist3(x, y, z, nbins)

if nargin > 3
  nbins = nbins;
else
  nbins = 35;
end

xbins = linspace(min(x), max(x), nbins);
ybins = linspace(min(y), max(y), nbins);
zbins = linspace(min(z), max(z), nbins);
D = zeros(nbins, nbins, nbins);

for i = 1:numel(x)
    xi = find((x(i) > xbins), 1, 'last');
    yi = find((y(i) > ybins), 1, 'last');
    zi = find((z(i) > zbins), 1, 'last');
    D(xi, yi, zi) = D(xi, yi, zi) + 1;
end

% [xn, xedge] = histcounts(x, nbins); 
% [yn, yedge] = histcounts(y, nbins); 
% [zn, zedge] = histcounts(z, nbins); 

figure;

%isosurface(D); colormap default; colorbar;
cmap = parula(256);
orthosliceViewer(D, 'Colormap', cmap);
set(gcf, 'Position',  [300, 300, 800, 800])
end
