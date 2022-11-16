fnames = [...
    'gemstrigger_2010072201.fid/fid';  % none
    'gemstrigger_2010072216.fid/fid';  % 1 hz
    'gemstrigger_2010072202.fid/fid';  % 10 hz
    'gemstrigger_2010072203.fid/fid';  % 500hz
    'gemstrigger_2010072204.fid/fid';  % 1000hz
    'gemstrigger_2010072205.fid/fid';  % 10,000 hz
    'gemstrigger_2010072206.fid/fid';  % 500, 000 hz
    'gemstrigger_2010072207.fid/fid';  % 1,000, 000 hz
    ];

ims = [];
for g = 1:size(fnames,1)
    
    filename = fnames(g,:);
    [kspace,np,nt] = ReadVarian2D(filename);
    kdat = fliplr(flipud(reshape(kspace,np/2,nt)));
    
%     framesize = size(kdat,2);
%     th=2*pi*[1:framesize]/4;
%     for ky=1:size(kdat,1)
%         kdat(ky,:)= kdat(ky,:) .* exp(-i*th);
%     end
    
    image_data = ifftshift(ifft2(fftshift(kdat)));
    tmp = image_data;
%     tmp(1:2:end, 2:2:end) = -tmp(1:2:end, 2:2:end);
    ims = [ims;  tmp(:)'];
end

g=1
while g < size(fnames,1) && g>0
    subplot(211)
    imagesc(reshape(abs(ims(g,:)), size(image_data))),colormap(gray);
    
    subplot(212)
    imagesc(reshape(angle(ims(g,:)), size(image_data))),colormap(gray);
    
    disp('use the LR arroow keys to move back and forth')
    [a b c] = ginput(1);
    switch(c)
        case 29
            g=g+1
        case 28
            g = g-1
    end
    
end