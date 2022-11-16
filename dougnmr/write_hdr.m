function write_hdr(name,hdr)

% function write_hdr(name,hdr)
% Luis hernandez
% last edit 1-7-98
%
% Writes the analyze format header file from a file 'name' 

% The function opens a file and writes the structure to it
% hdr = struct(...
%      'sizeof_hdr', fread(pFile, 1,'int32'),...
%      'pad1', setstr(fread(pFile, 28, 'char')),...
%      'extents', fread(pFile, 1,'int32'),...
%      'pad2', setstr(fread(pFile, 2, 'char')),...
%      'regular',setstr(fread(pFile, 1,'char')), ...
%      'pad3', setstr(fread(pFile,1, 'char')),...
%      'dims', fread(pFile, 1,'int16'),...
%      'xdim', fread(pFile, 1,'int16'),...
%      'ydim', fread(pFile, 1,'int16'),...
%      'zdim', fread(pFile, 1,'int16'),...
%      'tdim', fread(pFile, 1,'int16'),...
%      'pad4', setstr(fread(pFile,20, 'char')),...
%      'datatype', fread(pFile, 1,'int16'),...
%      'bits', fread(pFile, 1,'int16'),...
%      'pad5', setstr(fread(pFile, 6, 'char')),...
%      'xsize', fread(pFile, 1,'float'),...
%      'ysize', fread(pFile, 1,'float'),...
%      'zsize', fread(pFile, 1,'float'),...
%      'glmax', fread(pFile, 1,'int32'),...
%      'glmin', fread(pFile, 1,'int32'),... 
%      'descrip', setstr(fread(pFile, 80,'char')),...
%	'aux_file'        , setstr(fread(pFile,24,'char'))',...
%	'orient'          , fread(pFile,1,'char'),...
%				0 = transverse,1 = coronal, 2=sagittal
%	'origin'          , fread(pFile,5,'int16'),...
%	'generated'       , setstr(fread(pFile,10,'char'))',...
%	'scannum'         , setstr(fread(pFile,10,'char'))',...
%	'patient_id'      , setstr(fread(pFile,10,'char'))',...
%	'exp_date'        , setstr(fread(pFile,10,'char'))',...
%	'exp_time'        , setstr(fread(pFile,10,'char'))',...
%	'hist_un0'        , setstr(fread(pFile,3,'char'))',...
%	'views'           , fread(pFile,1,'int32'),...
%	'vols_added'      , fread(pFile,1,'int32'),...
%	'start_field'     , fread(pFile,1,'int32'),...
%	'field_skip'      , fread(pFile,1,'int32'),...
%	'omax'            , fread(pFile,1,'int32'),...
%	'omin'            , fread(pFile,1,'int32'),...
%	'smax'            , fread(pFile,1,'int32'),...
%	'smin'            , fread(pFile,1,'int32') );
%      )
   
%global SPM_scale_factor
SPM_scale_factor = 1;
   
   % Read in Headerfile into the hdrstruct
   [pFile,messg] = fopen(name, 'wb');
   if pFile == -1
      errormesg(messg);   
   end
   
        
      fwrite(pFile, hdr.sizeof_hdr,'int32');
      fwrite(pFile,hdr.pad1, 'char');
      fwrite(pFile,hdr.extents, 'int32');
      fwrite(pFile,hdr.pad2, 'char');
      fwrite(pFile,hdr.regular','char');
      fwrite(pFile,hdr.pad3', 'char');
      fwrite(pFile,hdr.dims','int16');
      fwrite(pFile,hdr.xdim','int16');
      fwrite(pFile,hdr.ydim', 'int16');
      fwrite(pFile,hdr.zdim', 'int16');
      fwrite(pFile,hdr.tdim', 'int16');
      fwrite(pFile,hdr.pad4', 'char');
      fwrite(pFile,hdr.datatype,'int16');
      fwrite(pFile,hdr.bits','int16');
      fwrite(pFile,hdr.pad5','char');
      fwrite(pFile,hdr.xsize', 'float');
      fwrite(pFile,hdr.ysize', 'float');
      fwrite(pFile,hdr.zsize', 'float');
      fwrite(pFile,hdr.pad6','char');
      fwrite(pFile,hdr.glmax', 'int32');
      fwrite(pFile,hdr.glmin', 'int32');
      fwrite(pFile,hdr.descrip','char');
      fwrite(pFile,hdr.pad6','char');
%      	fwrite(pFile,hdr.aux_file','char');
%      	fwrite(pFile,hdr.orient','char');
%	fwrite(pFile,hdr.origin','int16');      
%	fwrite(pFile,hdr.generated','char');
%	fwrite(pFile,hdr.scannum','char');
%	fwrite(pFile,hdr.patient_id','char');
%	fwrite(pFile,hdr.exp_date','char');
%	fwrite(pFile,hdr.exp_time','char');
%	fwrite(pFile,hdr.hist_un0','char');
%	fwrite(pFile,hdr.views', 'int32');
%	fwrite(pFile,hdr.vols_added', 'int32');
%	fwrite(pFile,hdr.start_field', 'int32');
%	fwrite(pFile,hdr.field_skip', 'int32');
%	fwrite(pFile,hdr.omax', 'int32');
%	fwrite(pFile,hdr.omin', 'int32');
%	fwrite(pFile,hdr.smax', 'int32');
%	fwrite(pFile,hdr.smin', 'int32');
%	

    fseek(pFile, 112, 'bof');
    fwrite(pFile, SPM_scale_factor , 'float');
      

   fclose(pFile);
   
return


