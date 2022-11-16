function hdr = read_hdr(name)
%
% (c) 2006 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
% function hdr = read_hdr(name)
% Luis hernandez
% last edit 8-24-2006
%
% Loads the analyze format header file from a file 'name' 
%
% The function returns a structure defined as 
% hdr = struct(...
% 0      'sizeof_hdr'	, fread(pFile, 1,'int32'),...
% 4      'pad1'		, setstr(fread(pFile, 28, 'char')),...
% 32     'extents'	, fread(pFile, 1,'int32'),...
% 36     'pad2'		, setstr(fread(pFile, 2, 'char')),...
% 38     'regular'	,setstr(fread(pFile, 1,'char')), ...
% 39     'pad3'		, setstr(fread(pFile,1, 'char')),...
% 40     'dims'		, fread(pFile, 1,'int16'),...
% 42     'xdim'		, fread(pFile, 1,'int16'),...
% 44     'ydim'		, fread(pFile, 1,'int16'),...
% 46     'zdim'		, fread(pFile, 1,'int16'),...
% 48     'tdim'		, fread(pFile, 1,'int16'),...
% 50     'pad4'		, setstr(fread(pFile,20, 'char')),...
% 70     'datatype'	, fread(pFile, 1,'int16'),...
% 72     'bits'		, fread(pFile, 1,'int16'),...
% 74     'pad5'		, setstr(fread(pFile, 6, 'char')),...
% 80     'xsize'		, fread(pFile, 1,'float'),...
% 84     'ysize'		, fread(pFile, 1,'float'),...
% 88     'zsize'		, fread(pFile, 1,'float'),...
% 92     'pad6'		, setstr(fread(pFile, 48, 'char'))...
% 140    'glmax'		, fread(pFile, 1,'int32'),...
% 144    'glmin'		, fread(pFile, 1,'int32'),... 
% 148    'descrip'	, setstr(fread(pFile, 80,'char')),...
% 228	   'aux_file'     , setstr(fread(pFile,24,'char'))',...
% 252    'orient'       , fread(pFile,1,'char'),...
% 253    'origin'       , fread(pFile,5,'int16'),...
% 263    'generated'    , setstr(fread(pFile,10,'char'))',...
% 273    'scannum'      , setstr(fread(pFile,10,'char'))',...
% 283    'patient_id'   , setstr(fread(pFile,10,'char'))',...
% 293    'exp_date'        , setstr(fread(pFile,10,'char'))',...
% 303    'exp_time'        , setstr(fread(pFile,10,'char'))',...
% 313    'hist_un0'        , setstr(fread(pFile,3,'char'))',...
% 316    'views'           , fread(pFile,1,'int32'),...
% 320    'vols_added'      , fread(pFile,1,'int32'),...
% 324    'start_field'     , fread(pFile,1,'int32'),...
% 328    'field_skip'      , fread(pFile,1,'int32'),...
% 332    'omax'            , fread(pFile,1,'int32'),...
% 336    'omin'            , fread(pFile,1,'int32'),...
% 340    'smax'            , fread(pFile,1,'int32'),...
% 344    'smin'            , fread(pFile,1,'int32') );
%

% $Id: read_hdr.m 1279 2014-03-24 20:06:25Z klitinas $
% $HeadURL: svn+ssh://klitinas@anger.engin.umich.edu/svn/matlab/img/trunk/read_hdr.m $

global SPM_scale_factor endian

isNIFTI = 0;

% checking to see that it's a .hdr file
suffix = name(end-3:end);
switch suffix
	case '.img'
		name = sprintf('%s.hdr',name(1:end-4));
	case '.hdr'
		name = sprintf('%s.hdr',name(1:end-4));
	% Next two cases accounts for NIFTI files
	case '.nii'
		isNIFTI = 1;
	case 'i.gz'
		fprintf('\nUnzipping ... %s', name);
		str = sprintf('!gunzip %s', name);
		eval(str);
    		name = name(1:end-3);
	otherwise
		name = sprintf('%s.hdr',name);
end

   % first detect which endian file we're opening
   % the first byte should be the number 348 (the size of the file)
   [pFile,messg] = fopen(name, 'r','native');
   if pFile == -1
      fprintf('Error opening header file: %s',name); 
      return;
   end

   tmp = fread(pFile,1,'int32');
   if strcmp(computer,'GLNX86') | ...
           strcmp(computer , 'PCWIN') | ...
           strcmp(computer,'GLNXA64') | ...
           strcmp(computer,'MACI') | ...
           strcmp(computer,'MACI64')
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
 
   % Now Read in Headerfile into the hdrstruct
   [pFile,messg] = fopen(name, 'r', endian);
   if pFile == -1
      sprintf('Error opening header file: %s - %s ',name, messg); 
      return;
   end

   % last check to see if we're dealing with a NIFTI file
   fseek(pFile, 344, 'bof');
   nifti_magic = char(fread(pFile, 4, 'char'));

   fseek(pFile,0,'bof'); 

   nstr = nifti_magic(1:3);
   if strcmp(nstr','ni1') | strcmp(nstr','n+1')
      isNIFTI=1;
   end

   if isNIFTI
       fclose(pFile);
       SPM_scale_factor =1;
       fprintf('\nFound %s : reading NIFTI  header %s\n', nifti_magic, name);
       h = read_nii_hdr(name);
       hdr = nii2avw_hdr(h);
       return
   end

   hdr = struct(...
      'sizeof_hdr'      , fread(pFile, 1,'int32'),...
      'pad1'            , setstr(fread(pFile, 28, 'char')),...
      'extents'         , fread(pFile, 1,'int32'),...
      'pad2'            , setstr(fread(pFile, 2, 'char')),...
      'regular'         , setstr(fread(pFile, 1,'char')), ...
      'pad3'            , setstr(fread(pFile,1, 'char')),...
      'dims'            , fread(pFile, 1,'int16'),...
      'xdim'            , fread(pFile, 1,'int16'),...
      'ydim'            , fread(pFile, 1,'int16'),...
      'zdim'            , fread(pFile, 1,'int16'),...
      'tdim'            , fread(pFile, 1,'int16'),...
      'pad4'            , setstr(fread(pFile,20, 'char')),...
      'datatype'        , fread(pFile, 1,'int16'),...
      'bits'            , fread(pFile, 1,'int16'),...
      'pad5'            , setstr(fread(pFile, 6, 'char')),...
      'xsize'           , fread(pFile, 1,'float'),...
      'ysize'           , fread(pFile, 1,'float'),...
      'zsize'           , fread(pFile, 1,'float'),...
      'pad6'            , setstr(fread(pFile, 48, 'char')),...
      'glmax'           , fread(pFile, 1,'int32'),...
      'glmin'           , fread(pFile, 1,'int32'),... 
      'descrip'         , setstr(fread(pFile, 80,'char')),...
      'aux_file'        , setstr(fread(pFile,24,'char'))',...
      'orient'          , fread(pFile,1 ,'char'),...
      'origin'          , fread(pFile,5,'int16'),...
      'generated'       , setstr(fread(pFile,10,'char'))',...
      'scannum'         , setstr(fread(pFile,10,'char'))',...
      'patient_id'      , setstr(fread(pFile,10,'char'))',...
      'exp_date'        , setstr(fread(pFile,10,'char'))',...
      'exp_time'        , setstr(fread(pFile,10,'char'))',...
      'hist_un0'        , setstr(fread(pFile,3,'char'))',...
      'views'           , fread(pFile,1,'int32'),...
      'vols_added'      , fread(pFile,1,'int32'),...
      'start_field'     , fread(pFile,1,'int32'),...
      'field_skip'      , fread(pFile,1,'int32'),...
      'omax'            , fread(pFile,1,'int32'),...
      'omin'            , fread(pFile,1,'int32'),...
      'smax'            , fread(pFile,1,'int32'),...
      'smin'            , fread(pFile,1,'int32') ... 
      );
  
    % this is where SPM hides its scaling factor when it writes floating point images in 
    % byte format to save space.  Those crafty guys ...
    fseek(pFile, 112, 'bof');
    SPM_scale_factor = fread(pFile, 1, 'float');
   

   fclose(pFile);

if strcmp(suffix, 'i.gz')
    fprintf('\nRe-zipping ... %s', name);
    str = sprintf('!gzip %s', name);
    eval(str);
end
        
return


