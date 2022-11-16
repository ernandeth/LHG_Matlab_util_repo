% process the data from the T1 and T2 measurements of Blood at 7T 
l = dir('t1_prof-03*.fid')
FID_filename = l.name;

[np,ns,nv,pss,lro,thk,data,seqfil,array, array_val]=ReadFID5(FID_filename);