function [array, hdr] = read_img2(hdr, name)

% [image_data, hdr] = read_img2([hdr,] name)
%
% Luis hernandez
% last edit 6-29-2000
%
% Loads the data from an analyze format file 'name' containing mutislice image data
% 		hdr: is a header structure (use read_hdr())
% 		name:  file name (.img)
% The function returns a three dimensional array 
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%

global SPM_scale_factor endian

if nargin==1
    % if there is only one argument, that means that it is
    % the file name
    name = hdr;    
end

if length(name)> 3
   suffix = name(end-3:end);
   if strcmp(suffix,'.img')
       hdr = read_hdr(sprintf('%s.hdr',name(1:end-4)));
   elseif strcmp(suffix,'.hdr')
       hdr = read_hdr(name);
       name = sprintf('%s.img',name(1:end-4));
   else
       hdr = read_hdr(sprintf('%s.hdr',name));
       name = sprintf('%s.img',name);
   end
else
    hdr = read_hdr(sprintf('%s.hdr',name));
    name = sprintf('%s.img',name);
end


if (abs(SPM_scale_factor < 0.000000000001)) 
    SPM_scale_factor=1;
end
%SPM_scale_factor


[pFile,messg] = fopen(name, 'r',endian);
if pFile == -1
    disp(messg);   
    return;
end
   

% figure out the format of the data file from the header.
% only three types are supported right now ...
switch hdr.datatype     
   case 0
      fmt = 'int8';
   case 2
      fmt = 'uint8';
   case 4
      fmt = 'short';
   case 8
      fmt = 'int';
   case 16
      fmt = 'float';
   case 32
      fmt = 'float';
      xdim = hdr.xdim * 2;
      ydim = hdr.ydim * 2;
   case 64
        fmt = 'int64';

   otherwise
         errormesg(sprintf('Data Type %d Unsupported. Aborting',hdr.datatype));
         return
end

	   
if hdr.tdim>1
    fprintf('Warning:  this file has a time series. Reading the first one ONLY');
end

% Create a 3D array with the data ...
array = zeros(hdr.xdim, hdr.ydim, hdr.zdim);
for i=1:hdr.zdim
	a  = (fread(pFile,[hdr.xdim, hdr.ydim], fmt)) ;
	array (:,:,i) = a;	
end

fclose(pFile);
array = array * SPM_scale_factor;

return

