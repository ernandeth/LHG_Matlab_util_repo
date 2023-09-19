function wtextfile(filename,var,fmt)
fid = fopen(filename,"w");
fprintf(fid, [fmt '\n'], var);
fclose(fid)