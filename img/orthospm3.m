function orthospm3(roi,threshold, onsets, window ,spm_file, anat_file, tseries_path, func_root )
% function result = orthospm3(  roi (Num. Neighbors),
%                               threshold , 
%                               onsets, 
%                               window 
%                               [,spm_file , 
%                                anat_file, 
%                                tseries_path, 
%                                func_root_name] )
%
% this function takes a file in analyze format and 
% displays orthogonal sections thru the planes interscting
% at x,y,z.  It then overlays a Statistical Parameter Map 
% thresholded at "threshold".
%
% additionally, this function allows you to extract a time series
% from a chose pixel, just by clicking on it.  The time series is saved
% in the file "tdata.dat" if you use the RIGHT MOUSE BUTTON.
% This file gets overwritten everytime you select a new pixel.
% 
% *** This version of the program calls event_avg.m 
% in order to average all the events into one.  No deconvolution!!
%
% onsets and window must be in scan units
%
% assumption: time series and SPM map have the same dimensions
%
% June 2005:  included an image of all the events so you can see the trends
% over the different trials.  No upsampling of data by default.
%
global SPM_scale_factor wscale w2scale
global SPM_scale_factor doGFilter doDetrend
global myfig
global ACTIVE
ACTIVE=1;
Evfig=figure;
set(gcf,'Position',[1 1 360,320]);
myfig=figure;

if nargin < 4
    [name path] = uigetfile('*.img','Select SPM *.img file');
    spm_file = strcat(path,name);
end

% figure out name for spm file
sz = size(spm_file);
imgname = strcat(spm_file(1,1:sz(2)-4) , '.img');
hdrname = strcat(spm_file(1,1:sz(2)-4) , '.hdr');
%label=sprintf('FUNCTIONAL:  ....%s',imgname(end-15:end));

spm_hdr = read_hdr(hdrname);
spm_data = read_img2(spm_hdr,imgname);
spm_data(find(spm_data==NaN))=0;
spm_scale =  SPM_scale_factor;

if nargin < 4
    [name path] = uigetfile('*.img','Select anatomical *.img file');
    anat_file= strcat(path,name)
end

% figure out name for anatomical file
sz = size(anat_file);

imgname = strcat(anat_file(1,1:sz(2)-4) , '.img');
hdrname = strcat(anat_file(1,1:sz(2)-4) , '.hdr');
%label=sprintf('%s\nANATOMY: ....%s', label, imgname(end-15:end));

hdr = read_hdr(hdrname);
anat_data = read_img2(hdr,imgname);

% interpolate the paramter map to fit the anatomical one
% note the transpose....
[x,y,z] = meshgrid(1:spm_hdr.ydim , 1:spm_hdr.xdim, 1:spm_hdr.zdim);
[xi,yi, zi] = meshgrid(1:hdr.ydim , 1:hdr.xdim, 1:hdr.zdim);

yi = yi * spm_hdr.ydim/hdr.ydim;
xi = xi * spm_hdr.xdim/hdr.xdim;
zi = zi * spm_hdr.zdim/hdr.zdim;

if spm_hdr.zdim > 1
    spm_data2 = interp3(x,y,z, spm_data, xi,yi,zi,'nearest');
else
    spm_data2=spm_data;
end
spm_data2(find(isnan(spm_data2)))=0;
%spm_data = spm_data2;

%whos

% threshold the map at the 50%
if nargin ==0
    threshold = mean(mean(mean(spm_data2))) * 3
end
spm_data2(find(spm_data2 <= threshold )) = 0;

% scale the maps to use the whole colormap
anat_data = anat_data * 256 / max(max(max(anat_data)));
spm_data2 = spm_data2 * 256 / max(max(max(spm_data2)));
    
if ~isempty(w2scale)
    spm_data2 = (spm_data2)*256 / w2scale(2) ;
else
    spm_data2 = spm_data2 * 256 / max(max(max(spm_data2)));
end

out_data = anat_data;
% manual scale
if ~isempty(wscale)
    out_data = (anat_data-wscale(1))*256 / wscale(2);
end
out_data(find(spm_data2)) =  256 + spm_data2(find(spm_data2))+1;


%configure the colormap:
set_func_colors


d=out_data;

%fprintf('\n first display the orthogonal sections');
x=ceil(hdr.xdim/2);
y=ceil(hdr.ydim/2);
z=ceil(hdr.zdim/2);

%colordef black 	
stretch = hdr.zsize/hdr.xsize;

%colordef black
[fig1, fig2, fig3] =  ov(hdr,d,x,y,z,0);

wholedata=[];
if nargin < 4
    % determine which files contain the time series..:
    f='dummmyname';
    file=[];
    path = [];
    op=pwd;
    while f ~= 0
        [f p] = uigetfile('*.img','Select a file from the next run.  cancel when finished ');
        
        if f~=0
            cd(p);
            
            fprintf('\nLoading the time series from ...%s / %s\n',p,  f(1:end-8))
            wholedata = [wholedata; read_img_series(f(1:end-8))];
            whos wholedata
            fprintf('\n...done\n')    
        end
    end
    cd(op);
    numRuns = size(f,1);
else
    op=pwd;
    cd(tseries_path);
    fprintf('\nLoading the time series from ...%s / %s\n',tseries_path,  func_root)
    wholedata = [wholedata; read_img_series(func_root)];
    
end

%tdata = xtract(spm_hdr, wholedata, [x-roi x+roi],[y-roi y+roi],[z-roi z+roi]);

[i j button] = ginput(1);
i=round(i);j=round(j);


%i=1;
while ACTIVE==1
    if i< -1, ACTIVE=0, end;
    fig = floor(gca);
    switch(fig)
        case floor(fig1)
            x=j;
            y=i;
        case floor(fig2)
            z=j;
            x=i;
        case floor(fig3)
            y=i;
            z=j;
    end
    

        
    xs = ceil(x*spm_hdr.xdim/hdr.xdim);
    ys = ceil(y*spm_hdr.ydim/hdr.ydim);
    zs = ceil(z*spm_hdr.zdim/hdr.zdim);
    
    str=sprintf('(x,y,z)=  (%3.2f %3.2f %3.2f) mm, \n (%3d %3d %3d)vx , \n val= %d',...
        hdr.xsize*(xs-spm_hdr.origin(1)), ...
        hdr.ysize*(ys-spm_hdr.origin(2)), ...
        hdr.zsize*(zs-spm_hdr.origin(3)), ...
        xs,ys,zs, ...
        spm_data(xs,ys,zs)*spm_scale )
    
    %colordef black
    fprintf('\n picked x y z value :')
    disp([ xs ys zs spm_data(xs,ys,zs)*spm_scale ])
         
    [fig1, fig2, fig3] =  ov(hdr,d,x,y,z,0);
    %subplot(221), title(fig1,str) %, xlabel(label)
    fprintf('\n%s' , str);
    
    tmp = xtract(spm_hdr, wholedata, [x-roi x+roi],[y-roi y+roi],[z-roi z+roi]);
    
    if doGFilter
        g = make_gaussian(50,4,100);
        
        disp('Gaussian filter ...')
        MeanBefore=mean(tmp);
        tmp = conv(g,tmp);
        tmp = tmp(50:end-50);
        MeanAfter=mean(tmp);
    end
    

    series_mean = mean(tmp);
    
    if doDetrend
        disp('detrending ...')
        dtmp = mydetrend(tmp');
        load coeffs.mat
        tmp=dtmp;
        tmp = 100 * tmp / series_mean;
    else
        tmp = 100 * (tmp' - series_mean) / series_mean;
    end
    %tdata = [tdata tmp];

    [ev_avg ev_std] = event_avg(tmp,onsets,window,1);
    
    
    
    %subplot 224, plot(events','y.');hold on
    %subplot 224, plot(10*spm_hrf(2),'y'),hold on
    %keyboard
    % subplot 224, plot(ev_avg,'r');axis tight ; hold off; 
    subplot 224, errorbar([1:length(ev_avg)], ev_avg,ev_std, ev_std);
    axis tight ; 
    hold off; 
    set(gca, 'Xtick',[0:2:window])
    
    events = [ev_avg' ev_std'];
    
    if button==3
        disp('saved avgResponse.dat and tdata.dat');
        save avgResponse.dat events -ASCII
        save tdata.dat tmp -ASCII
    end
    
    % now display the image of event averages    
    load allevents
    figure(Evfig);
    imagesc(allevents);colorbar
    xlabel('Time (scans)'), ylabel('Events')
    % return focus to other figure
    figure(myfig);
    
    [i j button] = ginput(1);
    i=round(i);j=round(j);
    
end 
colordef white
return
