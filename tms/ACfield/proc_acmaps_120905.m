%{  
% part 1: 

freq = 2e3;
GAMMA = 42.56 * 1e6 / 1e4;   % gyromagentic constant in Hz/Gauss
F = 6 * GAMMA/freq; 
% integrate 2 - cos(2*pi*f*t) from 0 to 3/f
threshold = 5000;

[data h] = fdf2nii('./gems_trigger2_20120905_04.img/',0);
imnums = [1:28];

% begin by reading in all the images
raw = [];
for n=1:14
	fname = sprintf('./gems_trigger2_20120905_%02d.fid/fid', imnums(n) )
	tmp = fidFT(fname);
	lightbox(abs(reshape(tmp,64,64,50)));
	raw = [raw; tmp(:)'];
	
end

% compute the Bz projections from the phase differences
% 
bz = angle( raw(2:2:end, :)  ./ raw(1:2:end,:) ) / F;

for pos=1:7
	% mask out the stuff outside the object
	msk = abs(raw(2*pos-1,:));
	msk(abs(msk) < threshold) = 0;
	msk(abs(msk) > 0) = 1;
	bz(pos,:) = msk.* bz(pos,:) / F; % note adjustmend into Gauss!

	% display the Bz maps to make sure:
	tmp = reshape(bz(pos,:),h.dim(2), h.dim(3), h.dim(4));
	lightbox(tmp,[-0.8 0.8],6);

    % write out the magnitude image and the Bz projection at each position
	% Note that the Bz fields are in milliGauss !!
	write_nii(sprintf('pos_%02d.nii', pos), abs(raw(pos*2,:)), h,0);
	write_nii(sprintf('Bz_%02d.nii', pos), 1000 * bz(pos,:), h,0);
	
end

% The next shell script corregisters all the images to the first one in the series
% it also produces the appropriate transformation matrix
!flirt_batch pos_01 pos_
!fslmerge -t r_allpositions rpos_*.nii

% applying the transformations to the phase maps
disp('Apply transformations to the phase images ...');
Bznames = dir('Bz_*');
MATnames = dir('./MATS/MAT*')


for f=1:7
   str = sprintf('!flirt -in %s -ref %s -applyxfm -init MATS/%s -out r%s', Bznames(f).name, Bznames(f).name, MATnames(f).name, Bznames(f).name)
   eval(str)
end;

!fslmerge -t r_allBz rBz_*.nii


% Make the projection matrix from the transformation matrices for each
% rotation.  First read the appropriate transformation matrices in
% and invert them to figure the right transformation to come back
cd MATS
mfiles = dir('MAT_*');
M = zeros(length(mfiles),3);
for ii=1:length(mfiles)
    tmp =load(mfiles(ii).name);
    itmp=inv(tmp);
    M(ii,:) = itmp(3,1:3);
end;
% Now split the matrix into two - recall that this data set was collected at two frequencies.
M2k = M;

iM2k=pinv(M2k);
disp(sprintf('condition number of aa is %f',cond(M2k)));

cd ..

% finally solve the problem.
% solve the inverse problem for Bx and By
% following the notation of Hernandez et al.

[tmp h] = read_nii_img('r_allBz.nii');

% read in the data.  They will be in miliGauss!!
Bz_measured = tmp;

Bx=zeros(1, size(Bz_measured,2));
By=zeros(1, size(Bz_measured,2));
Bz=zeros(1, size(Bz_measured,2));
dim = size(Bz_measured);

temp2=Bz_measured;
temp=iM2k*temp2;
Bx=temp(1,:);
By=temp(2,:);
Bz=temp(3,:);

B_mag=sqrt(Bx.^2+By.^2+Bz.^2);

h.dim(5) = 1;



% % mask out pixels that are too noisy
 msk = lightbox('pos_01.nii');
 msk(abs(msk) < threshold) = 0;
 msk(abs(msk) > 0) = 1;

Bxm=Bx.*msk(:)';
Bym=By.*msk(:)';
Bzm=Bz.*msk(:)';
B_magm=sqrt(Bxm.^2+Bym.^2+Bzm.^2);

write_nii('./Bx.nii',Bxm,h,0);
write_nii('./By.nii',Bym,h,0);
write_nii('./Bz.nii',Bzm,h,0);


% reshape the images so that I can look at them
Bxm = reshape(Bxm,h.dim(2), h.dim(3), h.dim(4));
Bym = reshape(Bym,h.dim(2), h.dim(3), h.dim(4));
Bzm = reshape(Bzm,h.dim(2), h.dim(3), h.dim(4));
B_magm = reshape(B_magm,h.dim(2), h.dim(3), h.dim(4));

save tms_wkspace 
%}
load tms_wkspace
B_magm = log(B_magm); 
fovx = h.dim(2)*h.pixdim(2);
fovz = h.dim(4)*h.pixdim(4);
fovy = fovx;

[xx yy zz]= meshgrid(linspace(0.05,fovx-0.05,64), linspace(0.05,fovy-0.05,64),linspace(0.05,fovz-0.05,50));
Bx = reshape(Bxm,64,64,50);
By = reshape(Bym,64,64,50);
Bz = reshape(Bzm,64,64,50);

% 
% [startx starty startz] = meshgrid(fovx*linspace(-0.2,0.2, 5), fovy*linspace(-0.2 ,0.2, 5), 0.25*fovz);
% 
% streamline(xx,yy,zz,Bx,By,Bz,startx,starty,startz);
% streamline(xx,yy,zz,-Bx,-By,-Bz,startx,starty,startz);
% 
% axis([-fovx/2 fovx/2 -fovy/2 fovy/2 -fovz/2 fovz/2])    


%figure
X = 32;Y=19; Z=38;

%for Y=22
	
myFOV = [X-18 X+18 Y-18 Y+18 Z-35 Z+5];


subplot(131)
cla
daspect([fovx fovy fovz])

streamslice(xx*64/fovx,yy*64/fovy,zz*50/fovz, Bx,By,Bz, [X], [], []);
hold on
slice(B_magm,X,1,1)
colormap(hot)
%caxis([ 0 12])
shading flat
view(90,0)
drawnow
line([X X],[Y Y],[1 50],'LineWidth',1,'LineStyle','--','Color','w')
line([X X],[1 64],[Z Z],'LineWidth',1,'LineStyle','--','Color','w')
axis(myFOV)
hold off
zlabel('Z'); 
ylabel('Y');
xlabel('X');

subplot(132)
cla
daspect([fovx fovy fovz]);
streamslice(xx*64/fovx, yy*64/fovy, zz*50/fovz, Bx,By,Bz, [], [Y], []);
hold on
slice(B_magm,1,Y,1)
colormap(hot)
%caxis([ 0 12])
shading flat
view(0,0)
line([X X],[Y Y],[1 50],'LineWidth',1,'LineStyle','--','Color','w')
line([1 64],[Y Y ],[Z Z],'LineWidth',1,'LineStyle','--','Color','w')
axis(myFOV)
drawnow
hold off
zlabel('Z'); 
ylabel('Y');
xlabel('X');


subplot(133)
cla
daspect([fovx fovy fovz]);
streamslice(xx*64/fovx,yy*64/fovy,zz*50/fovz, Bx,By,Bz, [], [], [Z]);
hold on
slice(B_magm,1,1,Z)
colormap(hot)
%caxis([ 0 12])
shading flat
view(0.01,90)
line([X X],[1 64],[Z Z],'LineWidth',1,'LineStyle','--','Color','w')
line([1 64],[Y Y ],[Z Z],'LineWidth',1,'LineStyle','--','Color','w')
axis(myFOV)
hold off
zlabel('Z'); 
ylabel('Y');
xlabel('X');


drawnow

%end

print -depsc ./stream_BxByBz

figure
subplot(131)
lightbox(Bx(:,:,22:41),[-10 10],[5]);, title('B_x');dofontsize(16);  colorbar off
subplot(132)
lightbox(By(:,:,22:41),[-10 10],[5]);, title('B_y');dofontsize(16);  colorbar off
subplot(133)
lightbox(Bz(:,:,22:41),[-10 10],[5]);, title('B_z');dofontsize(16); 
colormap hsv

%print -depsc ./lightbox_BxByBz




return




