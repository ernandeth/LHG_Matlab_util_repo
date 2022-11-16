function ts = timeseries2(cmd,imnames,O)
% function ts = timeseries2(cmd,imnames,option struct)
% *****************************************************
% cmd : 'slice' 'voxel' 'volume' 'multi'
% Imnames:  full path and filenames for each image.
%		one cell per session of images.  see spm_list_files or tor_list_files
%
% Options - enter in struct O.[field]
% O.include	=	number of images to include
%			should be cell array with vector of imgs to include for each session
% O.exclude	=	numbers of images to exclude in each session, cell array of length [# sessions]
% O.mask	=	masking image, use for volume method
% O.coord	=	voxel x,y,z coordinate to get, use for voxel method
% O.coords	=	voxel x,y,z coordinates to get, use for multi method
% O.slice	=	voxel z value of slice to get, use for slice method
% O.nohdrcheck	=	enter anything in this field to suppress header checking 
%			(otherwise checks for equal SPM scaling factors)
% O.verbose	=	enter anything to get verbose output.

% 'all' entered as last argument finds the ts for ALL the voxels in the region
% - normally finds the average of voxels in the region.
%
%
% voxel [i j k] is 3 column vectors. 
% i is rows of matrix, x coordinate in analyze img format
% j is cols of matrix, y coordinate in analyze
% k is slice number.
%
% ASSUMES BRAIN IS SIDEWAYS - so that x coord = col of matrix = y axis of actual brain.
% 'voxel' tested 1/29/01 for single_subj_T1.img.
% 'multi' tested same day, with s8t1.img (symmetrical x,y)
%
% Tor Wager, 10/17/01

format compact;
if ~iscell(imnames)
	tempNames{1} = imnames;
	imnames = tempNames;
end


% Set up image file related arguments
% ***************************************************************

file = 1;		% read from file rather than workspace; always true here.

% number of images and include images - shorten imnames if necessary
% ---------------------------------------------------------------
if isfield(O,'include'),
	if ~(length(imnames) == length(O.include)), 
		error('Img names and include vector must have equal numbers of cells.')
	end

	for i = 1:length(imnames)
		if any(O.include{i} > length(imnames{i})),
			error(['Session ' num2str(i) ': include image number exceeds number of images in session.'])
		end

		imnames{i} = imnames{i}(O.include{i},:);
	end

	
else 
	nimgs = 0;
	for i = 1:length(imnames)
		nimgs = nimgs + size(imnames{i},1);
	end
end



% exclude images if desired
% ---------------------------------------------------------------
if isfield(O,'exclude'),
	if ~(length(imnames) == length(O.exclude)), 
		error('Img names and exclude vector must have equal numbers of cells.')
	end

	for i = 1:length(imnames)
		if any(O.exclude{i} > length(imnames{i})),
			error(['Session ' num2str(i) ': exclude image number exceeds number of images in session.'])
		end

		try
			imnames{i}(O.exclude{i},:) = [];
		catch
			warning('Problem excluding images...?')
			imnames,O.exclude
		end
	end
end



 
% Read header files and check SPM scaling factors
% ---------------------------------------------------------------
firstim = deblank(imnames{1});
hdr_name = [firstim(1,(1:end-3)) 'hdr'];

try
	hdr = read_hdr(hdr_name);
catch
	disp('Can''t read header file! Returning header file name in ts')
	disp(hdr_name)
	disp(['Exist header name = ' num2str(exist(hdr_name))'])
	ts = hdr_name;
	return
end
	
if ~(isfield(O,'nohdrcheck'))
	% check to see if all SPM scaling factors are 1!!!  And list them in spmscale.
	index = 1;
	for i = 1:length(imnames)
		for j = 1:size(imnames{i},1)
			myimg = deblank(imnames{i}(j,:));
			hdr = read_hdr([myimg(1:end-3) 'hdr']);
			spmscale(index) = hdr.SPM_scale;
			index = index + 1;
		end
		if ~(mean(spmscale) == spmscale(1)),
			warning('Timeseries: not all spm scale factors are the same!!!'),
		end
   	end
end   

% concatenate image names
% ---------------------------------------------------------------
a = imnames{1};
for i = 2:length(imnames)
	a = str2mat(a,imnames{i});
end
imnames = a;
clear a;




% Set up masking and voxel-selection related arguments
% ***************************************************************


% threshold masking - not necessary?
% ---------------------------------------------------------------
if isfield(O,'mask'),
	mask = O.mask; domask = 1;
else
	mask = []; domask = 0;
end


% method specific arguments
% ---------------------------------------------------------------
switch cmd
	case 'slice'
		if isfield(O,'z'),
			z = O.z;
		else
			error('Must enter ''z'' field in option struct when using slice method.')
		end

        case 'voxel'
               	if isfield(O,'coord'),
			coord = O.coord;
		else
			error('Must enter ''coord'' field in option struct when using voxel method.')
		end
        case 'volume'
		if isfield(O,'mask'),
               		mask = double(O.mask);domask = 1;
		else
			error('Must enter ''mask'' field in option struct when using volume method.')
		end
        case 'multi'
		if isfield(O,'coords')
               		coords = O.coords;
		else
			error('Must enter ''coords'' field in option struct when using multi voxel method.')
		end
        otherwise error('unknown command in first argument: use ''slice'' ''voxel'' ''volume'' or ''multi''')
end  % end switch






% Set up dimensions and check data type
% ***************************************************************  
 
	% swap x and y dims!  xdim on image is ydim in matrix (rows)
	% now xdim, ydim, zdim are in terms of the image, which is sideways, not the brain.
	xdim = hdr.ydim;
 	ydim = hdr.xdim;
   	zdim = hdr.zdim;

	
   switch hdr.datatype     
   case 0
      fmt = 'int8';
	  bitsperbyte = 8;
   case 2
      fmt = 'uint8';
	  bitsperbyte = 8;
	  bytespervoxel = 1;
   case 4
      fmt = 'short';
	  bitsperbyte = 8;
	  bytespervoxel = 2;
   case 8
      fmt = 'int';
   case 16
      fmt = 'float';
   		bitsperbyte = 8;
	  bytespervoxel = 4;   
   case 32
      fmt = 'float';
      xdim = hdr.xdim * 2;
      ydim = hdr.ydim * 2
   otherwise
         error(sprintf('Data Type %d Unsupported. Aborting',hdr.datatype))
   end  % end switch


onevox = bytespervoxel * bitsperbyte;
clear ts;

if isfield(O,'verbose'),
   	disp(['Timeseries: format is ' num2str(fmt)])
   	disp(['bits per voxel = ' num2str(onevox) ', bytes per voxel = ' num2str(bytespervoxel) ', x,y dims = ' num2str(xdim) ' ' num2str(ydim)])
end







% Load the images and get the timeseries
% *************************************************************** 

switch cmd
   
case 'slice'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if nargin > 2
   	shift = bytespervoxel * (z-1) .* xdim .* ydim;
		%shift = onevox .* zshift;
   else error('must enter slice number as argument.'),end
   if file
      if isfield(O,'verbose'),disp([num2str(nimgs) ' files total']);,end
      t = cputime;
   	for i = 1:nimgs
   				fid = fopen(deblank(imnames(i,:)),'r');
   		try,
			status = fseek(fid,shift,-1);,
		catch,
			disp(['File is:' deblank(imnames(i,:))]),
			disp(['Number ' num2str(i) ' in series.'])
			warning(['file not opened properly.']),
			fid = fopen(deblank(imnames(i,:)),'r');
		end
         fclose(fid);
      
         if mod(i, 300) == 0 & isfield(O,'verbose')
            disp(['Done with file number ' num2str(i) '. ' num2str(cputime -t) 's since last update']);
            t = cputime;
         end
         	%if mod(i,50) == 0,disp(['done ' num2str(i) '...']),end
      end
  else 
      ts = wildcard(:,:,z,:);
  end 
  return
  
case 'voxel'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(O,'verbose')
	disp(['voxel = ' num2str(coord)])
	disp(['dimensions = ' num2str(xdim) ' ' num2str(ydim) ' ' num2str(zdim)])
end

	% Tested only with multi - not with single voxel.  Should work ok, though.
	% If analyze format image, columns of file are x in brain, rows are y in brain
	% Imaging a slice with readim2 transposes the rows and columns - so it looks like the brain is sideways
	% This may be due to how they are read in readim2.
	% Voxel coordinates are [x y z] in image, also [x y z] in brain if in Analyze format.
   	zshift = (coord(3)-1) .* xdim .* ydim;
   	yshift = (coord(2)-1) .* ydim;               % y should index # of rowsin file, cols in transposed image
   	xshift = (coord(1)-1);                       % x should index # of columns in file, rows in transposed image

if file
	shift = bytespervoxel * (zshift + yshift + xshift);	% shifts in bytes; 2 bytes per voxel
	for i = 1:nimgs
		fid = fopen(deblank(imnames(i,:)),'r');
   		try,
			status = fseek(fid,shift,-1);,
		catch,
			disp(['File is:' deblank(imnames(i,:))]),
			disp(['Number ' num2str(i) ' in series.'])
			warning(['file not opened properly.']),
			fid = fopen(deblank(imnames(i,:)),'r');
		end
		if i == 1 & isfield(O,'verbose'),
			disp(['reading bits ' num2str(ftell(fid)+1) ' - ' num2str(shift + onevox)]),
		end
		if ~(ftell(fid) == shift ),warning('coordinates exceeded end of file!'),end
   		ts(i,1) = fread(fid,1,fmt);
       fclose(fid);
       if mod(i, 300) == 0 & isfield(O,'verbose')
            disp(['Done with file number ' num2str(i) '. ' num2str(cputime -t) 's since last update']);
            t = cputime;
       end
end  
	
else
   error('Must enter file names in this version.')
end

case 'volume'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
if domask
   if isfield(O,'verbose'),disp('masking image with mask variable.'),end
   if nargin == 5
       % 'all' doesn't work - too much memory!  try uint8 or some scalefactor
       disp('function will return timeseries for EACH voxel.')   
       %pack
   else
       if isfield(O,'verbose'),disp('function will return timeseries with average of all voxels in mask.'),end
   end
   nvox = sum(sum(sum(mask)));
   for i = 1:size(mask,3)
      if sum(sum(mask(:,:,i))) > 0,maxslice = i;,end
   end
   for i = size(mask,3):-1:1;
      if sum(sum(mask(:,:,i))) > 0,minslice = i;,end
   end
   zshift = (minslice-1) .* xdim .* ydim;
   shift = 2 .* zshift;
else 
   if isfield(O,'verbose'),disp('no mask specified - finding all values'),end
   minslice = 1;
   maxslice = zdim;
   shift = 0;
end
zeroa = zeros(ydim,xdim);
for i = 1:nimgs
   if file
      		fid = fopen(deblank(imnames(i,:)),'r');
   		try,
			status = fseek(fid,shift,-1);,
		catch,
			disp(['File is:' deblank(imnames(i,:))]),
			disp(['Number ' num2str(i) ' in series.'])
			warning(['file not opened properly.']),
			fid = fopen(deblank(imnames(i,:)),'r');
		end
      for j=minslice:maxslice 
	     a  = (fread(fid,[ydim,xdim], fmt)) ;
         array (:,:,j) = a;
      end
      for j = maxslice+1:zdim
         array(:,:,j) = zeroa;
      end
      fclose(fid);
      if mod(i, 300) == 0 & isfield(O,'verbose')
            disp(['Done with file number ' num2str(i) '. ' num2str(cputime -t) 's since last update']);
            t = cputime;
      end
   else
      array = wildcard(:,:,:,i);
   end
   if domask,array = array .* mask;,end
   if nargin == 5,ts(:,:,:,i) = array;              % ALL ts 
   else ts(i,1) = sum(sum(sum(array))) / nvox;      % AVG ts
   end
   clear array;
end   

case 'multi'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
for k =  1:size(coords,1)
	% checked for accuracy 1/29/01, along with 'voxel'
    coord = coords(k,:);
	% Tested 10/28/01 with simulated data, Tor.  Works.
	% If analyze format image, columns of file are x in brain, rows are y in brain
	% Imaging a slice with readim2 transposes the rows and columns - so it looks like the brain is sideways
	% This may be due to how they are read in readim2.
   	zshift = (coord(3)-1) .* xdim .* ydim;
   	yshift = (coord(2)-1) .* ydim;               %%% y should index # of rws in file = cols in transposed image
   	xshift = (coord(1)-1);                       %%% x should index # of cols in file = rows in transposed image
    if file
	    shift = bytespervoxel .* (zshift + yshift + xshift);
	    for i = 1:nimgs
   		fid = fopen(deblank(imnames(i,:)),'r');
   		try,
			status = fseek(fid,shift,-1);,
		catch,
			disp(['File is:' deblank(imnames(i,:))]),
			disp(['Number ' num2str(i) ' in series.'])
			warning(['file not opened properly.']),
			fid = fopen(deblank(imnames(i,:)),'r');
		end
			if ~(ftell(fid) == shift ),warning('coordinates exceeded end of file!'),end
   		    ts.indiv(i,k) = fread(fid,1,fmt);
            fclose(fid);
	    if mod(i, 300) == 0 & isfield(O,'verbose')
            	disp(['Done with file number ' num2str(i) '. ' num2str(cputime -t) 's since last update']);
            	t = cputime;
            end
        end  
    else
        ts.indiv(:,k) = wildcard(coord(2),coord(1),coord(3),:);
    end
    if isfield(O,'verbose'),disp(['done ' num2str(k) ' voxels. coord = ' num2str(coords(k,:))]),end
end
ts.avg = mean(ts.indiv,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end					% end switch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return