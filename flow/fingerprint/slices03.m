function M=slices03(root, wscale)
M = [];
    [im h] = read_img(root);

for count=1:size(im,1)

    tmp=reshape(im(count,:) ,h.xdim,  h.ydim*h.zdim);
    %M = [M abs(tmp)'];
    M = [M (tmp)'];
end

if nargin==2
	show(M, wscale)
else
	show(M)
end

colormap gray(256)

xlabel('Time Points')
ylabel('slices')
