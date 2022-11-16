function result = read_cplx2(filename, linesize, nlines)
%function result = read_cplx2(filename, linesize, nlines );
% 
% This program reads in a file of complex point data, Line by line
% and returns a complex matrix.  The program always reads from the 
% end of the file, so you can skip the header ...
% ie - a set of spectra in a spectrospcopy experiment


result = [];
%[pf mesg]= fopen(filename,'r','ieee-be');
[pf mesg]= fopen(filename,'r','ieee-le');
if (pf == -1)
    errormesg(mesg)
end

fseek(pf,-linesize*nlines*2*4,'eof');


a = fread(pf, nlines*2*linesize,'long');
result = complex( a( 1:2:end) , a(2:2:end));
result = reshape(result, linesize, nlines);
fclose(pf)

return

