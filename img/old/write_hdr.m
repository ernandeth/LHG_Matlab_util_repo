function write_hdr(name,hdr)

% function write_hdr(name,hdr)
% Luis hernandez
% last edit 1-7-98
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%

% Writes the analyze format header file from a file 'name' 

% The function opens a file and writes the structure to it
% hdr = struct(...
%      'sizeof_hdr', fread(pFile, 1,'int32'),...
%      'pad1', setstr(fread(pFile, 28, 'uint8')),...
%      'extents', fread(pFile, 1,'int32'),...
%      'pad2', setstr(fread(pFile, 2, 'uint8')),...
%      'regular',setstr(fread(pFile, 1,'uint8')), ...
%      'pad3', setstr(fread(pFile,1, 'uint8')),...
%      'dims', fread(pFile, 1,'int16'),...
%      'xdim', fread(pFile, 1,'int16'),...
%      'ydim', fread(pFile, 1,'int16'),...
%      'zdim', fread(pFile, 1,'int16'),...
%      'tdim', fread(pFile, 1,'int16'),...
%      'pad4', setstr(fread(pFile,20, 'uint8')),...
%      'datatype', fread(pFile, 1,'int16'),...
%      'bits', fread(pFile, 1,'int16'),...
%      'pad5', setstr(fread(pFile, 6, 'uint8')),...
%      'xsize', fread(pFile, 1,'float'),...
%      'ysize', fread(pFile, 1,'float'),...
%      'zsize', fread(pFile, 1,'float'),...
%      'glmax', fread(pFile, 1,'int32'),...
%      'glmin', fread(pFile, 1,'int32'),... 
%      'descrip', setstr(fread(pFile, 80,'uint8')),...
%	'aux_file'        , setstr(fread(pFile,24,'uint8'))',...
%	'orient'          , fread(pFile,1,'uint8'),...
%				0 = transverse,1 = coronal, 2=sagittal
%	'origin'          , fread(pFile,5,'int16'),...
%	'generated'       , setstr(fread(pFile,10,'uint8'))',...
%	'scannum'         , setstr(fread(pFile,10,'uint8'))',...
%	'patient_id'      , setstr(fread(pFile,10,'uint8'))',...
%	'exp_date'        , setstr(fread(pFile,10,'uint8'))',...
%	'exp_time'        , setstr(fread(pFile,10,'uint8'))',...
%	'hist_un0'        , setstr(fread(pFile,3,'uint8'))',...
%	'views'           , fread(pFile,1,'int32'),...
%	'vols_added'      , fread(pFile,1,'int32'),...
%	'start_field'     , fread(pFile,1,'int32'),...
%	'field_skip'      , fread(pFile,1,'int32'),...
%	'omax'            , fread(pFile,1,'int32'),...
%	'omin'            , fread(pFile,1,'int32'),...
%	'smax'            , fread(pFile,1,'int32'),...
%	'smin'            , fread(pFile,1,'int32') );
%      )
   
global SPM_scale_factor
warning off
   
   % Read in Headerfile into the hdrstruct
   [pFile,messg] = fopen(name, 'wb');
   if pFile == -1
      errormesg(messg);   
   end
   
        
      fwrite(pFile, hdr.sizeof_hdr,'int32');
      fwrite(pFile,hdr.pad1, 'uint8');
      fwrite(pFile,hdr.extents, 'int32');
      fwrite(pFile,hdr.pad2, 'uint8');
      fwrite(pFile,hdr.regular','uint8');
      fwrite(pFile,hdr.pad3', 'uint8');
      fwrite(pFile,hdr.dims','int16');
      fwrite(pFile,hdr.xdim','int16');
      fwrite(pFile,hdr.ydim', 'int16');
      fwrite(pFile,hdr.zdim', 'int16');
      fwrite(pFile,hdr.tdim', 'int16');
      fwrite(pFile,hdr.pad4', 'uint8');
      fwrite(pFile,hdr.datatype,'int16');
      fwrite(pFile,hdr.bits','int16');
      fwrite(pFile,hdr.pad5','uint8');
      fwrite(pFile,hdr.xsize', 'float');
      fwrite(pFile,hdr.ysize', 'float');
      fwrite(pFile,hdr.zsize', 'float');
      fwrite(pFile,hdr.pad6','uint8');
      fwrite(pFile,hdr.glmax', 'int32');
      fwrite(pFile,hdr.glmin', 'int32');
      fwrite(pFile,hdr.descrip','uint8');
      fwrite(pFile,hdr.aux_file','uint8');
      fwrite(pFile,hdr.orient','uint8');
      fwrite(pFile,hdr.origin','int16');      
      fwrite(pFile,hdr.generated','uint8');
      fwrite(pFile,hdr.scannum','uint8');
      fwrite(pFile,hdr.patient_id','uint8');
      fwrite(pFile,hdr.exp_date','uint8');
      fwrite(pFile,hdr.exp_time','uint8');
      fwrite(pFile,hdr.hist_un0','uint8');
      fwrite(pFile,hdr.views', 'int32');
      fwrite(pFile,hdr.vols_added', 'int32');
      fwrite(pFile,hdr.start_field', 'int32');
      fwrite(pFile,hdr.field_skip', 'int32');
      fwrite(pFile,hdr.omax', 'int32');
      fwrite(pFile,hdr.omin', 'int32');
      fwrite(pFile,hdr.smax', 'int32');
      fwrite(pFile,hdr.smin', 'int32');
	

      fseek(pFile, 112, 'bof');
      fwrite(pFile, SPM_scale_factor , 'float');
      

      fclose(pFile);
   
return


