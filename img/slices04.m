function [M h]=slices04(root, wscale)
M = [];
    [im h] = readnii(root);

for count=1:size(im,4)

    tmp=reshape(im(:,:,:,count) ,h.dim(2),  h.dim(3)*h.dim(4));
    %M = [M abs(tmp)'];
    M = [M (tmp)'];
    
end

%M = abs(M);

if nargin==2
	imagesc(M)
    caxis(wscale)
else
	imagesc(M)
end


axis image
colormap gray(256)

xlabel('Time Points')
ylabel('slices')

colorbar