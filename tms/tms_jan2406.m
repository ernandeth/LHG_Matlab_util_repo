function tms_jan2406(datadir,P1,P2,P3,P4,P5,P6,xres)
%function tms_jan2406(datadir,P1,P2,P3,P4,P5,P6,xres)
% example : 
% tms_jan2406('/net/stanley/data/sangwool/tms/jan2406/','P26112.7','P26624.7','
% P27136.7','P27648.7','P28160.7','P28672.7',128);
%The following script is to reconstruct magnitude images and field maps using gsp21a. If any program fails, ask Luis about setting up paths properly. The necessary programs are : Luis' m file collection, SPM2 package, Sangwoo's m file collection
%addpath('~hernan/matlab/spm2b');
%addpath('~hernan/matlab/generic')

disp(sprintf('copying the P files in %s into current directory',datadir));
eval(['!cp ',datadir,P1,' Popen1']);
eval(['!cp ',datadir,P2,' Pon1']);
eval(['!cp ',datadir,P3,' Popen2']);
eval(['!cp ',datadir,P4,' Pon2']);
eval(['!cp ',datadir,P5,' Popen3']);
eval(['!cp ',datadir,P6,' Pon3']);
!rm -f *.img
!rm -f *.hdr
!rm -f ref*
!rm -f *.mat
!rm -f vol*
!rm -f rvol*
!rm -f sl*

%%% reconstruct the structural images and store them as analyze format in
%%% 'vol*.img' and 'vol*.hdr'
eval(sprintf('!/net/stanley/home/grlee/bin/gsp21a -n %d -fx -fy -m Popen3',xres));
eval(sprintf('!/net/stanley/home/grlee/bin/gsp21a -A -n %d -fx -fy -h Popen3',xres));
P = dir2names('vol_*.img');
aa=char(P(1));
eval(['!mv ',aa,' ',aa(1:22),'05',aa(25:end)]);
P = dir2names('vol_*.hdr');
aa=char(P(1));
eval(['!mv ',aa,' ',aa(1:22),'05',aa(25:end)]);

eval(sprintf('!/net/stanley/home/grlee/bin/gsp21a -n %d -fx -fy -m Popen2',xres)); 
eval(sprintf('!/net/stanley/home/grlee/bin/gsp21a -A -n %d -fx -fy -h Popen2',xres)); 
P = dir2names('vol_*.img');
aa=char(P(1));
eval(['!mv ',aa,' ',aa(1:22),'03',aa(25:end)]);
P = dir2names('vol_*.hdr');
aa=char(P(1));
eval(['!mv ',aa,' ',aa(1:22),'03',aa(25:end)]);

eval(sprintf('!/net/stanley/home/grlee/bin/gsp21a -n %d -fx -fy -m Popen1',xres));
eval(sprintf('!/net/stanley/home/grlee/bin/gsp21a -A -n %d -fx -fy -h Popen1',xres));
P = dir2names('vol_*.img');
aa=char(P(1));
eval(['!cp ',aa,' mask',aa(1:22),'01',aa(25:end)]);
P = dir2names('vol_*.hdr');
aa=char(P(1));
eval(['!cp ',aa,' mask',aa(1:22),'01',aa(25:end)]);

%%% reconstruct the open and on images and get the field maps and store
%%% them into an analyze format as 'open*.img' and 'open*.hdr'
mapdel=2500e-6;
eval(['!/net/stanley/home/grlee/bin/gsp21b -0 -p -fx -fy -n ' sprintf('%d',xres) ' Popen1']); 
P = dir2names('vol*.hdr');aa=char(P(1));
hdr1=read_hdr(aa);
hdr1.datatype=16;  % fake the header file
nsl = hdr1.zdim;
map=fmap(nsl,xres,mapdel);
write_img('open1.img',map,hdr1);
write_hdr('open1.hdr',hdr1);


mapdel=2500e-6;
eval(['!/net/stanley/home/grlee/bin/gsp21b -0 -p -fx -fy -n ' sprintf('%d',xres) ' Pon1']); 
map=fmap(nsl,xres,mapdel);
P = dir2names('vol*.hdr');aa=char(P(1));
hdr1=read_hdr(aa);
hdr1.datatype=16;  % fake the header file
write_img('on1.img',map,hdr1);
write_hdr('on1.hdr',hdr1);

mapdel=2500e-6;
eval(['!/net/stanley/home/grlee/bin/gsp21b -0 -p -fx -fy -n ' sprintf('%d',xres) ' Popen2']); 
map=fmap(nsl,xres,mapdel);
P = dir2names('vol*.hdr');aa=char(P(1));
hdr1=read_hdr(aa);
hdr1.datatype=16;  % fake the header file
write_img('open2.img',map,hdr1);
write_hdr('open2.hdr',hdr1);

eval(['!/net/stanley/home/grlee/bin/gsp21b -0 -p -fx -fy -n ' sprintf('%d',xres) ' Pon2']); 
map=fmap(nsl,xres,mapdel);
P = dir2names('vol*.hdr');aa=char(P(1));
hdr1=read_hdr(aa);
hdr1.datatype=16;  % fake the header file
write_img('on2.img',map,hdr1);
write_hdr('on2.hdr',hdr1);

eval(['!/net/stanley/home/grlee/bin/gsp21b -0 -p -fx -fy -n ' sprintf('%d',xres) ' Popen3']); 
map=fmap(nsl,xres,mapdel);
P = dir2names('vol*.hdr');aa=char(P(1));
hdr1=read_hdr(aa);
hdr1.datatype=16;  % fake the header file
write_img('open3.img',map,hdr1);
write_hdr('open3.hdr',hdr1);

eval(['!/net/stanley/home/grlee/bin/gsp21b -0 -p -fx -fy -n ' sprintf('%d',xres) ' Pon3']); 
map=fmap(nsl,xres,mapdel);
P = dir2names('vol*.hdr');aa=char(P(1));
hdr1=read_hdr(aa);
hdr1.datatype=16;  % fake the header file
write_img('on3.img',map,hdr1);
write_hdr('on3.hdr',hdr1);

!rm -f ref*
!rm -f sl*
%%
% at this stage, run SPM2 and realign magnitude images 
%starting with 'vol***.img'
P = dir2names('vol*.img');
spm_realign(P);
spm_reslice(P);

disp('run spm and check registration was done right.');
%keyboard;        
%% copy the realign matrices from structural (magnitude) images
Pr=dir2names('rvol*.mat');
P=dir2names('vol*.mat');
eval(['!cp ',char(Pr(1)),' open1.mat']);
eval(['!cp ',char(P(1)),' open2.mat']);
eval(['!cp ',char(P(2)),' open3.mat']);
eval(['!cp ',char(Pr(1)),' on1.mat']);
eval(['!cp ',char(P(1)),' on2.mat']);
eval(['!cp ',char(P(2)),' on3.mat']);

% mcflirt version:
! avwmerge -t tmp vol*.img
! mcflirt -in tmp -out rtmp -refvol 0 -cost normcorr -verbose 1 -stats -plots -mats        
! avwsplit rtmp.img
! rm rtmp.img rtmp.hdr tmp.hdr tmp.img
P = dir2names('rtmp.mat/');
eval(['!cp ',char(Pr(1)),' open1.mat']);
eval(['!cp ',char(P(1)),' open2.mat']);
eval(['!cp ',char(P(2)),' open3.mat']);
eval(['!cp ',char(Pr(1)),' on1.mat']);
eval(['!cp ',char(P(1)),' on2.mat']);
eval(['!cp ',char(P(2)),' on3.mat']);



% at this point, run 'SPM2 fmri' on matlab  and do the reslice only.
P = dir2names('open*.img');
spm_reslice(P);
P = dir2names('on*.img');
spm_reslice(P);

%% after running SPM2, (spm fmri)
%% read in the realigned field maps for on and open cases
ropen1=read_img2('ropen1.img');
ropen2=read_img2('ropen2.img');
ropen3=read_img2('ropen3.img');
ron1=read_img2('ron1.img');
ron2=read_img2('ron2.img');
ron3=read_img2('ron3.img');
% load in realigned-resliced magnitude images for mask
names = dir('rvol*.img');
str1=read_img2(names(1).name);
str2=read_img2(names(2).name);
str3=read_img2(names(3).name);
names=dir('rvol*.img');
maskim=read_img2(names(1).name);
%str11=phase_denoise_fun_3D(str1,ones(size(str1)),1);


% generate 3D mask
mask=maskim>max(maskim(:))*0.6;
%for ii=1:size(maskim,3)
%    mask(:,:,ii)=shrink(shrink(shrink(shrink(shrink(shrink(shrink(shrink(shrink(shrink(shrink(shrink(expand(expand(expand(expand(expand(expand(mask(:,:,ii)))))))))))))))))));
%end;
hdr1=read_hdr('open1.hdr');
write_img('bmask.img',mask,hdr1);
write_hdr('bmask.hdr',hdr1);
write_img('maskim.img',maskim,hdr1);
write_hdr('maskim.hdr',hdr1);

% calculate the induced Bz field from the phase difference map. 
mapdel=2.5e-3; % mapdelay in sec
gamma_bar=4258; % Hz / Gauss
Bz_1=[ron1-ropen1]/mapdel/gamma_bar*2*pi;
Bz_2=[ron2-ropen2]/mapdel/gamma_bar*2*pi;
Bz_3=[ron3-ropen3]/mapdel/gamma_bar*2*pi;

save Bzs Bz_1 Bz_2 Bz_3 mask
% read in the transformation matrix
load open1
M1=M(1:3,1:3);
load open2
M2=M(1:3,1:3);M2=inv(M2*inv(M1));
load open3
M3=M(1:3,1:3);M3=inv(M3*inv(M1));
M1=eye(3);

% finally solve the problem.
Bz_1(find(isnan(Bz_1)))=0;
Bz_2(find(isnan(Bz_2)))=0;
Bz_3(find(isnan(Bz_3)))=0;
% solve the inverse problem for Bx and By
% following the notation of Hernandez et al.
dim=size(Bz_1);
aa=inv([M1(3,:);M2(3,:);M3(3,:)]);
Bx=zeros(dim(1),dim(2),dim(3));
By=zeros(dim(1),dim(2),dim(3));
Bz=zeros(dim(1),dim(2),dim(3));
for ii=1:dim(1)
    for jj=1:dim(2)
        for kk=1:dim(3)
            if mask(ii,jj,kk)>=0
                temp=aa*[Bz_1(ii,jj,kk);Bz_2(ii,jj,kk);Bz_3(ii,jj,kk)];
                Bx(ii,jj,kk)=temp(1);By(ii,jj,kk)=temp(2);Bz(ii,jj,kk)=temp(3);
            end;
        end;
    end;
end;

%Bx=Bx.*mask;!cp Bz.mat Bx.mat
%!cp Bz.mat By.mat
%!cp Bz.img Bogus.img
%!cp Bz.hdr Bogus.hdr

%P = dir2names('B*.img');

%By=By.*mask;
%Bz=Bz.*mask;
B_mag=sqrt(Bx.^2+By.^2+Bz.^2);
if 0
for ii=170:-1:1
    subplot(231);quiver(sqz(Bx(:,:,ii)),sqz(Bz(:,:,ii)));title('xz');
    subplot(232);quiver(sqz(Bx(:,:,ii)),sqz(By(:,:,ii)));title('xy');
    subplot(233);quiver(sqz(By(:,:,ii)),sqz(Bz(:,:,ii)));title('yz');
    subplot(234);imagesc(flipud(sqz(B_mag(:,:,ii))),[0 50]);
    title(sprintf('magnitude Bz : axial slice number %d',ii));colorbar;
    subplot(236);imagesc(flipud(sqz(str1(:,:,ii))),[min(str1(:)), max(str1(:))]);
     title('magnitude image');pause;
end;
for ii=82:-1:1
    subplot(231);quiver(sqz(Bx(1:4:end,ii,1:4:end)),sqz(Bz(1:4:end,ii,1:4:end)));title('xz');
    subplot(232);quiver(sqz(Bx(1:4:end,ii,1:4:end)),sqz(By(1:4:end,ii,1:4:end)));title('xy');
    subplot(233);quiver(sqz(By(1:4:end,ii,1:4:end)),sqz(Bz(1:4:end,ii,1:4:end)));title('yz');
    subplot(234);imagesc(flipud(sqz(B_mag(1:1:end,ii,1:1:end))),[0 30]);
    title(sprintf('magnitude Bz : axial slice number %d',ii));colorbar;
    subplot(236);imagesc(flipud(sqz(str1(1:1:end,ii,1:1:end))),[min(str1(:)), max(str1(:))]);
     title('magnitude image');pause;
end;
for ii=85:-1:1
    subplot(231);quiver(sqz(Bx(ii,1:4:end,1:4:end)),sqz(Bz(ii,1:4:end,1:4:end)));title('xz');
    subplot(232);quiver(sqz(Bx(ii,1:4:end,1:4:end)),sqz(By(ii,1:4:end,1:4:end)));title('xy');
    subplot(233);quiver(sqz(By(ii,1:4:end,1:4:end)),sqz(Bz(ii,1:4:end,1:4:end)));title('yz');
    subplot(234);imagesc(flipud(sqz(B_mag(ii,1:1:end,1:1:end))),[0 30]);
    title(sprintf('magnitude Bz : axial slice number %d',ii));colorbar;
    subplot(236);imagesc(flipud(sqz(str1(ii,1:1:end,1:1:end))),[min(str1(:)), max(str1(:))]);
     title('magnitude image');pause;
end;
end;
hdr1=read_hdr('open1.hdr');
save results Bx By Bz B_mag mask

mask=B_mag<200;
%for ii=1:size(maskim,3)
%    mask(:,:,ii)=shrink(shrink(shrink(shrink(shrink(shrink(shrink(shrink(shrink(shrink(shrink(shrink(expand(expand(expand(expand(expand(expand(mask(:,:,ii)))))))))))))))))));
%end;
hdr1=read_hdr('open1.hdr');
write_img('bmask.img',mask,hdr1);
write_hdr('bmask.hdr',hdr1);
write_img('maskim.img',maskim,hdr1);
write_hdr('maskim.hdr',hdr1);

Bxm=Bx.*mask;
Bxm(abs(Bxm)>300)=0;
Bym=By.*mask;
Bym(abs(Bym)>300)=0;
Bzm=Bz.*mask;
Bzm(abs(Bzm)>300)=0;
hdr1=read_hdr('open1.hdr');

write_img('Bx.img',Bxm,hdr1);
write_img('By.img',Bym,hdr1);
write_img('Bz.img',Bzm,hdr1);
write_hdr('Bx.hdr',hdr1);
write_hdr('By.hdr',hdr1);
write_hdr('Bz.hdr',hdr1);
save results_masked_clipped Bxm Bym Bzm mask
keyboard;
%% after rotating Bx,By,Bz around manually using Bz in spm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!cp Bz.mat Bx.mat
!cp Bz.mat By.mat
!cp Bz.mat maskim.mat
!cp Bz.img Bogus.img
!cp Bz.hdr Bogus.hdr

P1 = dir2names('B*.img');
P2 = dir2names('maskim*.img');
P1(end+1)=P2;
spm_reslice(P1)



rBx=read_img2([],'rBx');
rBy=read_img2([],'rBy');
rBz=read_img2([],'rBz');
rmaskim=read_img2([],'rmaskim');

figure
subplot(311), hist(rBx(:),1000);
subplot(312), hist(rBy(:),1000);
subplot(313), hist(rBz(:),1000);

%%% run orthoq
orthoq('rB','rmaskim',[-30 30]);

mm=sqz(rmaskim(66,:,:));

for uu=69:69
imm(:,:,1)=fliplr(sqz(rBx(uu,:,:)).*mm);
imm(:,:,2)=fliplr(sqz(rBy(uu,:,:)).*mm);
imm(:,:,3)=fliplr(sqz(rBz(uu,:,:)).*mm);
im(imm,[-30 30]);
uu
pause;
mm=sqz(rmaskim(69,:,:)); mm=mm>max(mm(:))*0.1;
imm1=fliplr(sqz(rBz(69,:,:)).*mm);%coronal
mm=sqz(rmaskim(:,60,:)); mm=mm>max(mm(:))*0.1;
imm2=fliplr(sqz(rBz(:,60,:)).*mm);% sagittal
mm=sqz(rmaskim(:,:,9*16+3)); mm=mm>max(mm(:))*0.1;
imm3=fliplr(sqz(rBz(:,:,9*16+3)).*mm);%axial

h1=subplot(131);imagesc(imm1',[-30 30]);colormap gray;axis square;axis off;
hold on;plot(60,196-147,'y+');hold off;
set(h1,'Position',[ 0.1300   0.1100    0.2134    0.8150]);
h2=subplot(132);imagesc(imm2',[-30 30]);colormap gray;axis square;axis off;
hold on;plot(64-69,196-147,'y+');hold off;
set(h2,'Position',[ 0.1300+0.22   0.1100    0.2134    0.8150]);
h3=subplot(133);imagesc(imm3',[-30 30]);colormap gray;axis square;axis off;
hold on;plot(60,69,'y+');hold off;
set(h3,'Position',[ 0.1300+0.44   0.1100    0.2134    0.8150]);
colorbar
end;



