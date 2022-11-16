xres = 64;  slicenum=1; slice = randn(64,64); vnum = 1; Nslices=10; Nframes = 5;
volumes = zeros(Nframes, xres*xres*Nslices);

for count=1:Nslices*Nframes,
    slice = count* ones(64,64);
    begin = (slicenum -1)* xres *xres +1
    finish = (slicenum)* xres *xres
     volumes(vnum,  begin:finish ) = slice(:)' ;
 
     if (slicenum==Nslices), 
         lightbox(reshape(volumes(vnum,:),xres,xres,Nslices) );
         vnum = vnum+1;
         slicenum = 0;
         drawnow, 
     end;    
     slicenum = slicenum + 1
end