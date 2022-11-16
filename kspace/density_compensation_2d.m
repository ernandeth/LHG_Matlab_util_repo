function w = density_compensation_2d(k)
% function area = voronoidens(k);
%
% input:  k = kx + i ky is the  k-space trajectory
% output: area of cells for each point 
%           (if point doesn't have neighbors the area is NaN)
%
% Written by John Pauly, modified by Michael Lustig

% Thanks to this magical function, and help from Ray.
% Further modified for repeating k-points. 
% Shouchang Guo, 2017-08-24


r = max(abs(k(:)));
k = [k(:); r*1.005*exp(i*2*pi*[1:256]'/256)];
kx = real(k);
ky = imag(k);

% [row,column] = size(kx);
% uncomment these to plot voronoi diagram
% [vx, vy] = voronoi(kx,ky);
% figure, plot(kx,ky,'r.',vx,vy,'b-'); axis equal

fkxy = [kx(:),ky(:)];
[kxy,~,ik] = unique(fkxy, 'rows', 'stable'); %%

% returns vertices and cells of voronoi diagram
[V,C] = voronoin(kxy);
area = [];
for j = 1:length(kxy)
    if ~isempty(C{j})
        x = V(C{j},1); y = V(C{j},2); lxy = length(x);
        A = abs(sum(0.5*(x([2:lxy 1])-x(:)).*(y([2:lxy 1])+y(:))));
    else
        A = inf;
    end
    area = [area A];
end

ik = ik(1:end-256);
area = area(1:end-256); 

area = area(:)/sum(area(:));
w = area(ik);
[bincounts, ind] = histc(ik, unique(ik));
w = w./bincounts(ind);
w = w/max(abs(w(:)));




