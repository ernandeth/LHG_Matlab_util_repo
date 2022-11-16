function orthospm8(onsets, window ,anat_file, tseries_path, func_root,voxelFile )
% function result = orthospm8(
%                          onsets, 
%                          window 
%                          anat_file, 
%                          tseries_path, 
%                          func_root_name,
%                          voxelFile )
%
%  replaces(orthospm6)
%
% this function takes a file in analyze format and 
% displays orthogonal sections
%
% the data extracted are those from the voxels specified in the 
% TEXT file "voxelFile". 
%
% WARNING : the voxels must be in the same size as those from the 
% anat file for things to display correctly!
%
% ** this function is NOT interactive on purpose **
%
% It saves the files:
% tdata.data - the time series extracted
% avg.dat - the averaged events
% pixels.dat - the xyz coords of the pixels used 
%  
% *** This version of the program calls event_avg.m 
% in order to average all the events into one.  No deconvolution!!
%
% this version loads all the data into memory first.
%
% onsets and window must be in scan units
%
% The program replaces the first few points as disdaqs with the 
% value of the next point.  
% Change the value of disdaqs if you want to change that !!
%
% use the global "upsample" to determine interpolation in case 
% the onset times are non-integer

disdaqs = 10
global SPM_scale_factor upsample
%close all
figure

% figure out name for anatomical file
sz = size(anat_file);
imgname = strcat(anat_file(1,1:sz(2)-4) , '.img');
hdrname = strcat(anat_file(1,1:sz(2)-4) , '.hdr');
%label=sprintf('%s\nANATOMY: ....%s', label, imgname(end-15:end));

hdr = read_hdr(hdrname);
anat_data = read_img2(hdr,imgname);


cx = round(hdr.xdim/2);
cy = round(hdr.ydim/2);      
cz = round(hdr.zdim/2);      

colordef black
mygray = [0:1/127:1]' * [1 1 1]; 
colormap (mygray);

d=anat_data;
[fig1, fig2, fig3] =  ov(hdr,d,cx,cy,cz,0);

%if button==3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the list of voxels 

roi_xyz = load(voxelFile);          

% Make sure there are voxels above threshold
if size(roi_xyz) ~= [0 0]  
    
    % make a gaussian filter here	
    g = make_gaussian(50,4,100);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nLoading the time series ...%s\n', file(1:end-8))
    wholedata = read_img_series(file(1:end-8));
    whos wholedata
    fprintf('\n...done\n')    
         
    tdata=[];
    %tmp = xtractlist(hdr, wholedata, roi_xyz);
    
    fprintf('calculating the appropriate voxel coordinates...\n')
    inds = sub2ind([hdr.xdim hdr.ydim hdr.zdim],roi_xyz(:,1),roi_xyz(:,2), roi_xyz(:,3));
    num=0;
    
    tcourse=zeros(size(wholedata,1),1);
    
    for n=1:size(list,1)
        tcourse = tcourse + ...
            wholedata( :,inds(n) );
        num = num +1;
    end
    
    tmp = tcourse/num;

    if (disdaqs>0) tmp(1:disdaqs) = tmp(disdaqs);  % get rid of the pesky disdaqs!
          
        save raw.dat tmp -ASCII
        % filtering stuff:
        % tmp = smoothdata2(tmp, TR, 0.009, 0.3, 3); 
        tmp = mydetrend(tmp');
        % load the mean (DC term) and divide by it
        % (compute percentage)
        load coeffs.mat
        mm = coeffs(end);
        tmp = 100 * tmp / mm;
        
        tmp = conv(g,tmp);
        tmp = tmp(50:end-50);
        tdata = [tdata tmp];
        
        if isempty(upsample)
            [ev_avg ev_std] = event_avg(tdata,onsets,window,1);
        else
            [ev_avg ev_std] = event_avg(tdata,onsets,window,upsample);
        end
        subplot(223), hold on, plot(roi_xyz(:,2), roi_xyz(:,1),'g*');
        subplot(222), hold on, plot(roi_xyz(:,1), roi_xyz(:,3),'g*');
        subplot(221), hold on, plot(roi_xyz(:,2), roi_xyz(:,3),'g*');
        
        subplot 224, plot(ev_avg,'r');axis tight ; hold off; 
        subplot(224), title (sprintf('ROI size: %d vox',size(roi_xyz,1)))
        set(gca, 'Xtick',[0:2:window])
        
        
        tmp = [ev_avg ev_std];
        save avg.dat tmp -ASCII
        save tdata.dat tdata -ASCII
        save voxels.dat roi_xyz -ASCII
    end
    return
