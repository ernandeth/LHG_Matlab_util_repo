function swap_hdrx(name)
%
% function swap_hdrx(name)
%
% last edit 8.24.06  - Luis Hernandez-Garcia at UM
%
% this program will read in an analyze or nifti header
% it will then swap the xsize field (float at byte 80)
% for its absolute value.
% this is to counter act the effects of the FSL programs
% that stick a negative sign into the xsize of all the headers.
%
% this program does not touch anything else in the header!
%

   % first detect which endian file we're opening
   % the first byte should be the number 348 (the size of the file)
   [pFile,messg] = fopen(name, 'r','native');
   if pFile == -1
      msgbox(sprintf('Error opening header file: %s',name)); 
      return;
   end

   tmp = fread(pFile,1,'int32');
   if strcmp(computer,'GLNX86') | strcmp(computer , 'PCWIN')
       
       if tmp==348
           endian='ieee-le';
       else
           endian='ieee-be';
       end
       
   else
       if tmp==348
           endian='ieee-be';
       else
           endian='ieee-le';
       end
       
   end
   fclose(pFile);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   [pFile,messg] = fopen(name, 'r', endian);
   if pFile == -1
      sprintf('Error opening header file: %s - %s ',name, messg); 
      return;
   end
   hdr = fread(pFile,348,'uint8');
   fseek(pFile, 80, 'bof');
   xs = fread(pFile,1,'float');
   fclose(pFile);

   fprintf('\n%s Found %f : writing %f',name,  xs, abs(xs));
   [pFile,messg] = fopen(name, 'wb', endian);
   if pFile == -1
      sprintf('Error opening header file: %s - %s ',name, messg); 
      return;
   end
   fwrite(pFile,hdr,'uint8');
   fseek(pFile, 80, 'bof');
   fwrite(pFile, abs(xs), 'float');
   fclose(pFile);

return
        


