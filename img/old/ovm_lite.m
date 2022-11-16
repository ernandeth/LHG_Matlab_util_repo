function ovm_lite(threshold,spm_file, spm_file2, anat_file )
% function ovm_lite(threshold,spm_file, spm_file2, anat_file )
%
% this function overlays two activation maps on top of an anatomical 
% image and shows the overlaps in green at the specified threshold
% 
% the display is in orthogonal sections and it's NOT interactive
% (for use with batch jobs)
%
% The program writes a file called "overlaps.dat" with the xyz 
% voxel coordinates of the overlapping voxels
%
global SPM_scale_factor 
    sz = size(spm_file);
    
    imgname = strcat(spm_file(1,1:sz(2)-4) , '.img');
    hdrname = strcat(spm_file(1,1:sz(2)-4) , '.hdr');
    
    
    spm_hdr = read_hdr(hdrname);
    spm_data = read_img2(spm_hdr,imgname);
    spm_data(find(spm_data==NaN))=0;
    spm_scale =  SPM_scale_factor;

    sz = size(spm_file2);
    
    imgname = strcat(spm_file2(1,1:sz(2)-4) , '.img');
    hdrname = strcat(spm_file2(1,1:sz(2)-4) , '.hdr');
    
    
    spm_hdr2 = read_hdr(hdrname);
    spm_data2 = read_img2(spm_hdr2,imgname);
    spm_data2(find(spm_data2==NaN))=0;
    spm_scale2 =  SPM_scale_factor;

    voxList = find( (spm_data2 >= threshold) & (spm_data >=threshold));
    [vlx vly vlz] = ind2sub(size(spm_data), voxList);
    voxList = [vlx vly vlz]
    save overlaps.dat voxList -ASCII
    
    sz = size(anat_file);
    
    imgname = strcat(anat_file(1,1:sz(2)-4) , '.img');
    hdrname = strcat(anat_file(1,1:sz(2)-4) , '.hdr');
    
    hdr = read_hdr(hdrname);
    anat_data = read_img2(hdr,imgname);
    % interpolate the paramter map to fit the anatomical one
	% note the transpose....A weird quirk of meshgrid
    [x,y,z] = meshgrid(1:spm_hdr.ydim , 1:spm_hdr.xdim, 1:spm_hdr.zdim);
    [xi,yi, zi] = meshgrid(1:hdr.ydim , 1:hdr.xdim, 1:hdr.zdim);

    x = ( x - spm_hdr.origin(1) )* spm_hdr.xsize;
    y = ( y - spm_hdr.origin(2) )* spm_hdr.ysize;
    z = ( z - spm_hdr.origin(3) )* spm_hdr.zsize;
    xi = ( xi - hdr.origin(1) )* hdr.xsize;
    yi = ( yi - hdr.origin(2) )* hdr.ysize;
    zi = ( zi - hdr.origin(3) )* hdr.zsize;


    whos anat_data spm_data spm_data2
    hdr.origin'
    spm_hdr.origin'

    %whos x y z xi yi zi spm_data
    spm_dataR = interp3(x,y,z, spm_data, xi,yi,zi,'nearest');
    R = spm_dataR;
    spm_dataR(find(isnan(spm_dataR)))=0;
    %spm_data = spm_data2;
    spm_dataB = interp3(x,y,z, spm_data2, xi,yi,zi,'nearest');
    B = spm_dataB;
    spm_dataB(find(isnan(spm_dataB)))=0;
    
    
    % threshold the map at the 50%
    spm_dataR(find(spm_dataR <= threshold )) = 0;
    spm_dataB(find(spm_dataB <= threshold )) = 0;
    
    % scale the maps to use the whole colormap
    anat_data = anat_data * 2^7 / max(max(max(anat_data)));
    spm_dataR = spm_dataR * 2^7 / max(max(max(spm_dataR)));
    spm_dataB = spm_dataB * 2^7 / max(max(max(spm_dataB)));
    
    out_data = anat_data;
    out_data(find(spm_dataR)) =  2^7 + spm_dataR(find(spm_dataR))-1;
    out_data(find(spm_dataB)) =  2^8 + spm_dataB(find(spm_dataB))-1;
    out_data(find(spm_dataB & spm_dataR)) =  2^9 + spm_dataB(find(spm_dataB & spm_dataR))-1;

    %configure the colormap:
	figure
	mygray = [0:1/127:1]' * [1 1 1]; 
    
	myhot = [0:3:125]' * [1 0 0] ;
	tmp =   [0:3:125]' * [0 1 0] ;
	myhot = [myhot; tmp];
	tmp =   [0:3:125]' * [0 0 1];
	myhot =  [myhot;  tmp]/128;

	myhot(round(128/3): 128, 1) = 1;
	myhot(round(128*2/3):128,2) = 1;

	myblue = myhot;
	myblue(:,1) = myhot(:,3);
	myblue(:,3) = myhot(:,1);

	mygreen = myhot;
	mygreen(:,3) = 0;%mygray(:,1);
	mygreen(:,1) = 0;
	mygreen(:,2) = 1;

 	mymap = [mygray; myhot ; myblue ; mygreen];
    	colormap(mymap)


	axis off
    
        
    d=out_data;

  % display the orthogonal sections
    x = round(hdr.xdim/2);
    y = round(hdr.ydim/2);
    z = round (hdr.zdim/2);
    
    %colordef black 	
    stretch = hdr.zsize/hdr.xsize;
    
    xs=round(x);ys=round(y);zs=round(z);

	% calculate voxel coords. for spm data point
	xs = round((x - hdr.origin(1))*hdr.xsize/spm_hdr.xsize + spm_hdr.origin(1));
        ys = round((y - hdr.origin(2))*hdr.ysize/spm_hdr.ysize + spm_hdr.origin(2));
        zs = round((z - hdr.origin(3))*hdr.zsize/spm_hdr.zsize + spm_hdr.origin(3));

	% display the data:
	[fig1, fig2, fig3] =  ov(hdr,d,x,y,z,0);

	str=sprintf('(x,y,z)=  (%3.2f %3.2f %3.2f) mm, \n (x,y,z)=(%3d %3d %3d)vx\n R-val= %6.2f   B-val= %6.2f',...
		hdr.xsize*( x - hdr.origin(1)), ...
		hdr.ysize*( y - hdr.origin(2)), ...
		hdr.zsize*( z - hdr.origin(3)), ...
		xs,ys,zs, ...
		R(x,y,z)*spm_scale,...
		B(x,y,z)*spm_scale2 );
	subplot(221),title(str);      
	if isempty(voxList)
		subplot(223), title('No Overlaps!')
	end
		

     
    %colordef white
    return












