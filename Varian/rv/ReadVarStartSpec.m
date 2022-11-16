function[final_kspace,np,ntraces,nblocks]= ReadVarStartSpec(nspec)
%This file reads in Varian raw fid file 
%Output:final_kspace,np,ntraces,nblocks sw sfrq

%Example:
%[final_kspace,np,ntraces,nblocks]= ReadVarStart(3);

global pathname  np ntraces nblocks arraydim seqcon sw sfrq

if ischar(pathname)
    [filename,pathname]=uigetfile('*.*','Select fid file for spectroscopy data',pathname);
else
[filename,pathname]=uigetfile('*.*','Select fid file for the spectroscopy data');
end

if filename ~=0
    %Open the input file
    filename=[pathname filename];
    [fid,msg]=fopen(filename,'rt');
       if fid<0
        % There is an error
        str='File Open Failed'
       else
        filename
       end
else
    error('File Open Stopped')
    final_kspace=0;
    np=0;
    ntraces=0;
    nblocks=0;
end

 [kspace,np,ntraces,nblocks] = ReadVarian2D(filename);
 
%seqcon='nccnn' %this is for gems
%seqcon='ncsnn'; %this is for sems
filenameProc=[pathname 'procpar']
   [fid,msg]=fopen(filenameProc,'rt');
    if fid<0
        
     %There is an error
     str=['File Open Failed'];
     errordlg(str,title,'modal');
     
     else
   
     sw = ReadProcpar('sw',filenameProc)
     sfrq = ReadProcpar('sfrq',filenameProc)
     disp(sprintf('Reading in Varian Procpar files from %s',filename));
     arraydim=ReadProcpar('arraydim',filenameProc)
     seqcon = ReadProcpar('seqcon',filenameProc)
    end  
    
    
    final_kspace=kspace;
    
    if seqcon(2)&&seqcon(3)~='n'
        disp('this is csi data')
        %rearrange the data
        a=reshape(final_kspace,[nblocks np/2 ntraces]);
        z1=fftshift(fftn(a));
        figure;pcolor((squeeze(abs(z1(:,100,:)))));shading flat;colormap(gray)
    else
        
    figure;plot(real(final_kspace(nspec,:)));
    title('This is FID #1')
    figure;plot(real(fliplr(fftshift(fft(final_kspace(nspec,:))))));
    title('This is SPECTRUM #1')
        
    end
 
 
