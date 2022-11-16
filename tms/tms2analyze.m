function tms2analyze(tmsdata)
%
% function tms2analyze(tmsdata)
%
% Luis Hernandez
% last edit 8-12-98
%
% Converts TMS position and intensity data to Analyze format using
% a user selected analyze header file.
%
% generates an ANALYZE format slice set.
% Zeros are placed where no SPM values are speciefied.

   
   % Determine the dimensions of the spm->analyze output data
   % Read in Analyze hdr info:
   [f p] = uigetfile('*.hdr','Select Image file with appropriate dimensions - CANCEL to use defaults');
   hdrname = strcat(p,f);
   
   if hdrname == 0
      hdrname = '/export/home/hernan/matlab/default.hdr';
      errormesg('using the default header .../export/home/hernan/matlab/default.hdr');
   end

   % Extract the analyze header information

   hdr = read_hdr(hdrname);
      
   % Scale the xyz coordinates to pixel numbers 
   % create an array of zeros of the dimensions of the 
   % output file.
   % replace zeros at those locations with the parameter values
   % convert from cm to mm
   
   tmsdata(:,1) = tmsdata(:,1)*10 / hdr.xsize;
   tmsdata(:,2) = tmsdata(:,2)*10 / hdr.ysize;
   tmsdata(:,3) = tmsdata(:,3)*10 / hdr.zsize;
   
   img = zeros( 1, hdr.xdim* hdr.ydim* hdr.zdim);
   %size(img)

    
   pixels = ceil(tmsdata(:,1:3));
   sz = size(pixels)
   
   
   % Compute scaling factor for intensity datat
   if max(tmsdata(:,4)) ~= 0
      scale = 255 / max(tmsdata(:,4) );
   else
      scale = 0;
   end
   
   
   for i = 1 : sz(1)
      
     x = pixels(i,1);
	y = pixels(i,2);
	z = pixels(i,3);
     
     
     if (x>0) & (x <=hdr.xdim) & ...
          	(y>0) & (y <=hdr.ydim) & ...
		(y>0) & (z <=hdr.zdim)

		index = x + (y-1)* hdr.xdim + (z-1)*hdr.xdim*hdr.ydim;
	     
	     % Scale the TMS data so it can be displayed properly
	     if scale ~= 0
	        img(index) = tmsdata(i,4) * scale;
	     else
	        img(index) = 200;
	     end
          
       end
       
     
   end

   % Get the output filenames from the user
   [f p] = uiputfile('*.img','Save Analyze format file As ...');
   imgname = strcat(p,f);
   sz = size(imgname);
   hdrname = strcat(imgname(1:sz(2)-4) , '.hdr');
   
   % Write the header file
   disp(sprintf('writing %s',hdrname))
   
   hdr.bits = 8;
   hdr.datatype=2;
   
   write_hdr(hdrname,hdr);
   
   % Write the img file
   [fp mesg]= fopen(imgname,'wb');
   if fp ==-1
      errormesg(mesg);
      return;
   end
   
   
   switch hdr.bits      
	case 8
         fmt = 'uint8';
      case 16
         fmt = 'uint16';
      case 32
         fmt = 'uint32';
      otherwise
         errormesg(sprintf('Data Type %d Unsupported. Aborting',hdr.bits));
         return
      end

   fwrite(fp, img, fmt);
   disp(sprintf('writing %s',imgname));

   fclose(fp);
   
return
   
   













