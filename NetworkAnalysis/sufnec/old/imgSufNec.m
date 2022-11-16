function  imgSufNec(filename, xyz, DesMat)
%function imgSufNec(filename, [x,y,z], DesMat)
%
% (c) 2013 Luis Hernandez-Garcia
%      Alejandro Veloz-Magnus
% University of Michigan
% report bugs to:  hernan@umich.edu
%



[pth, name, suffix] = fileparts(filename);
[data, h] = read_img(filename);

% if suffix=='.nii'
% 	h = nii2avw_hdr(h);
% end 

% stuff to decorrelate from data
if nargin>2
    data = data - data*pinv(DesMat);
end

% Now binarize everything (detect events)
bdata = zeros(size(data));
for r=1:size(data,2)
    
    if data(1,r)~=0 && ~isnan(data(1,r))
        bdata(:,r) = detect_event(data(:,r));
        fprintf('\rDetecting Activation with MLE ...  %d', r);
    else
        bdata(:,r) = nan;
        
    end
end
write_img('Events.img', bdata*1e4, h);

% Extract the seed time course:
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

A = xtractor(h, bdata, x, y, z);

Suf = zeros(1,size(bdata,2));
Nec = zeros(1,size(bdata,2));

for r=1:size(bdata,2)
    B = bdata(:,r);
    
    if ~isnan(B(1))  && B(1)~=0
        [Nec(r), Suf(r)] = nec_suf(A, B);
    else
        Nec(r) = nan;
        Suf(r) = nan;
    end
    
end



h.tdim = 1;

global SPM_scale_factor 
SPM_scale_factor= 1e-4;

write_img('Sufficiency.img', Suf*1e4, h);
write_img('Necessity.img', Nec*1e4 , h);

return


function tcourse = xtractor(hdr, data, x, y, z)
%
% function tcourse = xtractROI(hdr, data, x, y, z)
%
% here the data is a 2D matrix with all the time points in rows
% x,y,z are vectors defining an arbitratry ROI
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
%

num = length(x);


tcourse=zeros(size(data,1),1);
inds = sub2ind([hdr.xdim, hdr.ydim, hdr.zdim], x , y, z);
raw = data(:,inds);
for count=1:length(tcourse)
    tmp = raw(count,:);
    tcourse(count) = mean(tmp(isfinite(tmp)));
end

fprintf('... %d voxels in ROI ', num)


return
