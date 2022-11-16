% Function that reads in binary .fid file created by Vnmrj.
%  [dta, fheader, bheader] = readfid(filename,n);
%    inputs:  filename: vnmr fid file to read
%             n: number of points at the end of fid to use for baseline correction (-1 = no correction)
%    Outputs: dta: returned fid
%             fheader: returned fid file header information
%             bheader: returned fid block header information for the last
%             block that was read.

function [dta,fheader,bheader]=readfid(filename,n)

   %filename: vnmr fid file to read
   %n: number of points at the end of fid to use for baseline correction
   %(-1 = no correction)
   %dta: returned fid
   %fheader: returned fid file header information
   %bheader: returned fid block header information for the last block that
   %was read

   fheader=struct('nblocks',{0},'ntraces',{0},'np',{0},'ebytes',{0},'tbytes',{0},'bbytes',{0},'vers_id',{0},'status',{0},'nb_headers',{0});
   bheader=struct('scale',{0},'status',{0},'index',{0},'mode',{0},'ctcount',{0},'lpval',{0},'rpval',{0},'lvl',{0},'tlt',{0});

   %Read vnmr file

   f=fopen(filename,'r','b');

   fheader.nblocks=fread(f,1,'long');  % 4 bytes, 32 bit
   fheader.ntraces=fread(f,1,'long');  % 4 bytes, 32 bit
   fheader.np=fread(f,1,'long');       % 4 bytes, 32 bit
   fheader.ebytes=fread(f,1,'long');   % 4 bytes, 32 bit
   fheader.tbytes=fread(f,1,'long');   % 4 bytes, 32 bit
   fheader.bbytes=fread(f,1,'long');   % 4 bytes, 32 bit
   fheader.vers_id=fread(f,1,'short'); % 2 bytes, 16 bit
   fheader.status=fread(f,1,'short');  % 2 bytes, 16 bit
   fheader.nb_headers=fread(f,1,'long'); % 4 bytes, 32 bit

   dta=[];

  for fv=1:fheader.nblocks
     bheader.scale=fread(f,1,'short');
     bheader.status=fread(f,1,'short');
     bheader.index=fread(f,1,'short');
     bheader.mode=fread(f,1,'short');
     bheader.ctcount=fread(f,1,'long');
     bheader.lpval=fread(f,1,'float');
     bheader.rpval=fread(f,1,'float');
     bheader.lvl=fread(f,1,'float');
     bheader.tlt=fread(f,1,'float');
    for nt = 1:fheader.ntraces
%     fprintf('Reading block %d\n',fv);
     d=fread(f,fheader.np,'float');

     d=transpose(reshape(d,2,size(d,1)/2));
     d=complex(d(:,1),d(:,2));
     if n>=0
        d=d-mean(d(end-n:end));
        dta=[dta,d];
     else
%          if nt == fheader.ntraces-4
%              dbstop readfid at 51
%          end
         dta=[dta,d];
     end
    end
     
  end
   
fclose(f);
end