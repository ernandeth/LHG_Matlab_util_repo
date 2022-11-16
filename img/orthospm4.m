function orthospm4(roi, threshold, onsets, window ,spm_file, anat_file, tseries_path, func_root )
% function result = orthospm4(
%                         roi_radius,
%                         threshold , 
%                         onsets, 
%                         window 
%                         [,spm_file , 
%                          anat_file, 
%                          tseries_path, 
%                          func_root_name] )
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
% The program replaces the first "disdaqs" (5) points with the 
% value of the next point.  
% Change the value of disdaqs if you want to change that !!
%

disdaqs = 5 
 
global SPM_scale_factor 
close all

    if nargin < 5
        [name path] = uigetfile('*.img','Select SPM *.img file');
        spm_file = strcat(path,name)
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

    if nargin < 5
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
    [x,y,z] = meshgrid(1:spm_hdr.ydim , 1:spm_hdr.xdim, 1:spm_hdr.zdim);
    [xi,yi, zi] = meshgrid(1:hdr.ydim , 1:hdr.xdim, 1:hdr.zdim);
    
    % optional in case the origin is not set in the image headers...
    %
    % spm_hdr.origin = [0 0 0];
    % hdr.origin = [0 0 0];
    %
    
    x = ( x - spm_hdr.origin(1) )* spm_hdr.xsize;
    y = ( y - spm_hdr.origin(2) )* spm_hdr.ysize;
    z = ( z - spm_hdr.origin(3) )* spm_hdr.zsize;
    
    xi = ( xi - hdr.origin(1) )* hdr.xsize;
    yi = ( yi - hdr.origin(2) )* hdr.ysize;
    zi = ( zi - hdr.origin(3) )* hdr.zsize;

    % yi = yi * spm_hdr.ydim/hdr.ydim;
    % xi = xi * spm_hdr.xdim/hdr.xdim;
    % zi = zi * spm_hdr.zdim/hdr.zdim;

    spm_data2 = interp3(x,y,z, spm_data, xi,yi,zi,'nearest');
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
    
    out_data = anat_data;
    out_data(find(spm_data2)) =  256 + spm_data2(find(spm_data2));
    
    %configure the colormap:
	set_func_colors

    
        
    d=out_data;

   fprintf('\nfirst display the orthogonal sections');
    x=ceil(hdr.xdim/2);
    y=ceil(hdr.ydim/2);
    z=ceil(hdr.zdim/2);
    
    %colordef black 	
    stretch = hdr.zsize/hdr.xsize;
    
    %colordef black
	[fig1, fig2, fig3] =  ov(hdr,d,x,y,z,0);
    
      
    if nargin < 5
        % determine which files contain the time series..:
        f='dummmyname';
        file=[];
        path = [];
        op=pwd;
        while f ~= 0
            [f p] = uigetfile('*.img','Select a file from the next run.  cancel when finished ');
            if f~=0
                file = [file ; f];
                path = [path ; p];
                cd(p);
            end
        end
        cd(op);
        numRuns = size(file,1);
    end
    
    [i j button] = ginput(1);
    i=round(i);j=round(j);
    
   
    %i=1;
    while i >= -10
        
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
        
        
        str=sprintf('(x,y,z)=  (%3.2f %3.2f %3.2f) mm,  \n (%3d %3d %3d)vx , \n val= %d',...
        	hdr.xsize*(xs-spm_hdr.origin(1)), ...
        	hdr.ysize*(ys-spm_hdr.origin(2)), ...
        	hdr.zsize*(zs-spm_hdr.origin(3)), ...
		    xs,ys,zs, ...
            spm_data(xs,ys,zs)*spm_scale );
        
        %colordef black
        [ x y z] 
        [fig1, fig2, fig3] =  ov(hdr,d,x,y,z,0);
        subplot(221), title(str) %, xlabel(label)
        
        

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find a list of voxels within the sphere VOI
        
        nlist = (2*roi+1)^3; % total # of voxels within a cube
        [xx,yy,zz] = ndgrid([xs-roi:xs+roi], [ys-roi:ys+roi], [zs-roi:zs+roi]); 
        xx_1d = reshape(xx,nlist,1); 
        yy_1d = reshape(yy,nlist,1); 
        zz_1d = reshape(zz,nlist,1); 
        
        cube_xyz = [xx_1d yy_1d zz_1d];   % a list of all voxels within the cube
        dist_cube = sqrt(...
            (hdr.xsize * (cube_xyz(:,1) - xs)).^2 + ...
            (hdr.ysize * (cube_xyz(:,2) - ys)).^2 + ...
            (hdr.zsize * (cube_xyz(:,3) - zs)).^2);    % a list of distances in mm
        
        sphere_xyz = cube_xyz(dist_cube <=roi, :); % a list of voxels within sphere
        num_voxels = size(sphere_xyz,1); % len_list is the # of voxels within sphere
        
        % list of values in the spm_data 
        roi_vals= [];
        for v=1: num_voxels
            roi_vals = [roi_vals; 
                spm_data2(sphere_xyz(v,1), sphere_xyz(v,2), sphere_xyz(v,3)) ];
        end
        roi_xyz = sphere_xyz( roi_vals >= threshold , :);          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make sure there are voxels above threshold
        if size(roi_xyz) ~= [0 0]  
           
            if nargin > 5
                numRuns=1;
                path=tseries_path;
                file=func_root;
            end
            % make a gaussian filter here	
            g = make_gaussian(50,4,100);
	           
            tdata=[];
            for n=1:numRuns
                tmp = timeplot4(path(n,:), file(n,:), roi_xyz);
                if (disdaqs>0) tmp(1:disdaqs) = tmp(disdaqs);  % get rid of the pesky disdaqs!
                save sprintf('raw_%d.dat',n) tmp -ASCII 

			
				% filtering stuff:
				% tmp = smoothdata2(tmp, TR, 0.009, 0.3, 3); 
				tmp = mydetrend(tmp');
				tmp = conv(g,tmp');
				tmp = tmp(50:end-50);
		
                tdata = [tdata tmp];
            end
            
            [ev_avg ev_std] = event_avg(tdata,onsets,window,10);
            
            subplot(221), hold on, plot(roi_xyz(:,2), roi_xyz(:,1),'g*');
            subplot(222), hold on, plot(roi_xyz(:,1), roi_xyz(:,3),'g*');
            subplot(223), hold on, plot(roi_xyz(:,2), roi_xyz(:,3),'g*');
           
            subplot 224, plot(ev_avg,'r');axis tight ; hold off; 
            subplot(224), title (sprintf('ROI at %d %d %d (%d vox)', xs , ys, zs, size(roi_xyz(2))))
            set(gca, 'Xtick',[0:2:window])
            
            if button==3
                tmp = [ev_avg ev_std];
                save avg.dat tmp -ASCII
                save tdata.dat tdata -ASCII
                save voxels.dat roi_xyz -ASCII
            end
        else
            subplot(224), title (sprintf('No active voxels at %d %d %d', xs , ys, zs))
        end
        [i j button] = ginput(1);
        i=round(i);j=round(j);
        
    end 
    %colordef white
    return
