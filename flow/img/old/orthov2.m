function result = orthov2( roi,root, tsroot, x, y, z)
% function result = orthov2( roi_size, root,   x, y, z)
% ...or    result=orthov2()
%
% The ROI is a cube defined by the specified number of voxels
% (roi_size) in each direction
%
% this function thakes a file in analyze format and 
% displays orthogonal sections thru the planes interscting
% at x,y,z
%
% additionally, this functiona allows you to extract a time series
% from a chose pixel, just by clicking on it.  The time series is saved
% in the file "tdata.dat" if you use the RIGHT MOUSE BUTTON.
% This file gets overwritten everytime you select a new pixel.
%
global fig
fig=figure;
global wscale

switch nargin
    case 0         % GUI selection of analtomical and time series
        roi=0;
        
        % select the image to display
        [name path] = uigetfile('*.img','Select Anatomical *.img file');
        name = strcat(path,name)
        
        sz = size(name);
        
        imgname = strcat(name(1,1:sz(2)-4) , '.img');
        hdrname = strcat(name(1,1:sz(2)-4) , '.hdr');
        h = read_hdr(hdrname);
        d = read_img2(h, imgname);
        
        % now select the time series
        [file path] = uigetfile('*.img','Select one of the Analyze files in the series ');
        
    case 1         % GUI selection of analtomical and time series
        
        % select the image to display
        [name path] = uigetfile('*.img','Select Anatomical *.img file');
        name = strcat(path,name)
        
        sz = size(name);
        
        imgname = strcat(name(1,1:sz(2)-4) , '.img');
        hdrname = strcat(name(1,1:sz(2)-4) , '.hdr');
        h = read_hdr(hdrname);
        d = read_img2(h, imgname);
        
        % now select the time series
        [file path] = uigetfile('*.img','Select one of the Analyze files in the series ');
        
    case 2  % NO GUI selection : time series and anatomical are the same.
        if root(end-3:end)=='.hdr'
            root=root(1:end-4);
        elseif root(end-3:end)=='.img'    
            root=root(1:end-4);
        end
        % image to display ....
        h = read_hdr(strcat(root,'.hdr'));
        d = read_img2(h, (strcat(root,'.img') ));
        
        % time series
        file = sprintf('%s.img',root);
        
    case 3  % NO GUI selection : time series and anatomical are the same.
        if root(end-3:end)=='.hdr'
            root=root(1:end-4);
        elseif root(end-3:end)=='.img'    
            root=root(1:end-4);
        end
        % image to display ....
        h = read_hdr(strcat(root,'.hdr'));
        d = read_img2(h, (strcat(root,'.img') ));
        
        % time series
        if tsroot(end-3:end)=='.hdr'
            tsroot=tsroot(1:end-4);
        elseif tsroot(end-3:end)=='.img'    
            tsroot=tsroot(1:end-4);
        end
        file = sprintf('%s.img',tsroot);
        
        
    otherwise
        if root(end-3:end)=='.hdr'
            root=root(1:end-4);
        elseif root(end-3:end)=='.img'    
            root=root(1:end-4);
        end
        % image to display ....
        h = read_hdr(strcat(root,'.hdr'));
        d = read_img2(h, (strcat(root,'.img') ));
        
        % time series
        file = sprintf('%s.img',root);
        
end

if nargin <5
    x=ceil(h.xdim/2);
    y=ceil(h.ydim/2);
    z=ceil(h.zdim/2);
end


stretch = h.zsize/h.xsize;

% display the orthogonal sections
colordef black 	
colormap(gray)

%     fprintf('\n(x,y,z)=  %d %d %d , val= %d, roi radius=%d', x, y, z, d(x,y,z), roi);
%     
%     colordef black
%     fig1=subplot (221), imagesc(squeeze(d(:,:,z))), axis ([1 h.ydim 1 h.xdim]) ,axis xy 
%     hold on; plot(y,x,'go');hold off;
%     fig2=subplot (222), imagesc(squeeze(d(:,y,:))'), axis ([1 h.xdim 0 h.zdim*stretch]),axis xy 
%     hold on; plot(x,z,'go');hold off;
%     fig3=subplot (223), imagesc(squeeze(d(x,:,:))'), axis ([1 h.ydim 0 h.zdim*stretch]), axis xy
%     hold on; plot(y,z,'go');hold off;

colordef black


%%%%%%
fprintf('\nLoading the time series ...%s\n', file(1:end-8))
wholedata = read_img_series(file(1:end-8));
whos wholedata
fprintf('\n...done\n')    
tdata = xtract(h, wholedata, [x-roi x+roi],[y-roi y+roi],[z-roi z+roi]);
%%%%%%%%%%%%%%%%%%%%%%%



subplot 224, plot(tdata)%, axis tight;

%%%%%%

k=0;

my_map=(0:255)';
my_map=[my_map my_map my_map]/256;
colormap(my_map);
dd = d*256/max(max(max(d)));

% scale image to fit colormap
range= max(max(max(d))) - min(min(min(d)));
dd = (d-min(min(min(d))))*256/range;
if ~isempty(wscale)
    dd = (d-wscale(1))*256 / wscale(2);
end

if ~isempty(wscale)
    dd = (d-wscale(1))*256 / wscale(2);
end


while k >= -10
    
    [fig1, fig2, fig3] =  ov(h,dd,x,y,z,roi);
    str = sprintf('\n(x,y,z)=  (%d %d %d) , val= %6.2f  \n', x, y, z, d(x,y,z));
    fprintf('\n%s' , str);
    
    tdata = xtract(h, wholedata, [x-roi x+roi],[y-roi y+roi],[z-roi z+roi]);
    
    
     subplot 224, plot(tdata)%, axis tight;
%     grid on
    
    [k j button] = ginput(1);
    k=round(k);j=round(j);
    fig = floor(gca);
    switch(fig)
        case floor(fig1)
            x=j;
            y=k;
        case floor(fig2)
            z=j;
            x=k;
        case floor(fig3)
            y=k;
            z=j;
    end
    
    if button==3
        save tdata.dat tdata -ASCII
    end
    % exiting the program when you click outside bounds
    if (k<=-1)   
        colordef white
        return
    end
end 
result = tdata;

return
