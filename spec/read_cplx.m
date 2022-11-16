function result = read_cplx(filename, linesize, nlines,hdrskip)
%function result = read_cplx(filename, linesize, nlines [,hdrskip]);
% 
% This program reads in a file of complex point data, Line by line
% and returns a complex matrix.  The program always reads the 
% end of the file, so you can skip the header ...
% ie - a set of spectra in a spectrospcopy experiment
%  
% Note that spectra headers are usually 48176 bytes

if nargin==4
    skip=hdrskip;
else
    skip=0;
end
result = [];
% [pf mesg]= fopen(filename,'r','ieee-be');
[pf mesg]= fopen(filename,'r','ieee-le');
if (pf == -1)
    errormesg(mesg)
end

if nargin==4
    fseek(pf,skip,'bof');
else
   fseek(pf,-linesize*nlines*2*2,'eof');
end


for i=1:nlines
	a = fread(pf,[2 linesize],'long');  % this is LONG
	b = complex(a(1, :) , a(2,:) )';
	result = [result  b];
end
fclose(pf)

return

