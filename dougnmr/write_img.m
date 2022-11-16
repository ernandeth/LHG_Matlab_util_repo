function  write_img(name, data, hdr)

% Luis hernandez
% last edit 4-24-2004
%
% function  write_img_data(name, data, hdr)
%
% Writes the data to an analyze format file 'name' containing mutislice image data 
% this also handles a timeseries in one file (ie -each image is a row )
%



   [pFile,messg] = fopen(name, 'wb');
   if pFile == -1
      errormesg(messg);   
      return;
   end
   
      
   switch hdr.datatype     
   case 2
      fmt = 'char';
   case 4
      fmt = 'short';
   case 8
      fmt = 'int';
   case 16
      fmt = 'float';
   case 32
      fmt = 'float';
           
   otherwise
      errormesg(sprintf('Data Type %d Unsupported. Aborting',hdr.datatype));
      return
      
   end

   
   if hdr.tdim>1
	fwrite(pFile, data', fmt);
	else
	fwrite(pFile, data, fmt);
   end

   fclose(pFile);
   
   
 return

