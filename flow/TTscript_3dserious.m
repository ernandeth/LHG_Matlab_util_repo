
Pfiles=dir('P*')

all = [];
%
for count=1:1 %length(Pfiles)
    
    % do the recon
    spiral3drecon_201119(Pfiles(count).name, 1, 0, 1);
    
    hdr = read_nii_hdr('timeseries_mag.nii');
    aslsub( 'timeseries_mag' , 1, 3, hdr.dim(5), 0,0,0);
    [s h]= read_img('mean_sub');
    all = [all ; s(:)'];
    str = ['!mv timeseries_mag.nii img_' Pfiles(count).name '.nii'];
    eval(str)

end
%}
%%
Pfiles=dir('img_P*');

all = [];
allcon = [];
alltag = [];
for count=1:length(Pfiles)
    
    fname = Pfiles(count).name;
    hdr = read_nii_hdr(fname);
    aslsub( fname(1:end-4) , 1, 5, hdr.dim(5), 0,1,0);
    
    [s h]= read_img('mean_sub');
    all = [all ; s(:)'];
    
    [c h]= read_img('mean_con'); 
    allcon = [allcon ; c(:)'];
    
    [t h]= read_img('mean_tag');
    alltag = [alltag ; t(:)'];
    
    str = ['!mv timeseries_mag.nii img_' Pfiles(count).name '.nii'];
    eval(str)

end
h.tdim = length(Pfiles);

write_img('all.img', all, h);
write_img('all_con.img', allcon, h);
write_img('all_tag.img', alltag, h);

close all

figure
slices02('all',[-250 250]) ;
axis image

hdr = read_hdr('all.hdr');
axis([1 hdr.xdim*hdr.tdim (hdr.zdim/2-2)*hdr.ydim (hdr.zdim/2+2)*hdr.ydim]);

figure
slices02('all_con') ;
axis image

hdr = read_hdr('all_con.hdr');
axis([1 hdr.xdim*hdr.tdim (hdr.zdim/2-2)*hdr.ydim (hdr.zdim/2+2)*hdr.ydim]);


figure
slices02('all_tag') ;
axis image

hdr = read_hdr('all_con.hdr');
axis([1 hdr.xdim*hdr.tdim (hdr.zdim/2-2)*hdr.ydim (hdr.zdim/2+2)*hdr.ydim]);
