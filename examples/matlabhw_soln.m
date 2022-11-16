% This assignment consists of a tutorial that will expose you to the 
% key concepts in matlab that you typically use in FMRI time series
% analysis
% 
% this assignment needs the following files in the working directory:
% read_hdr.m, read_img2.m, spm_hrf.m, tmp.hdr, tmp.img, SPM_fMRIDesMtx.mat and two
% directories full of images:  SPM, timeseries.  Please verufy that you
% have those befpre beginging the exercise.
%

% Topic 1:  Variables and their types
% Type these in line by line and paste the output below the line. Note how matlab
% behaves. You will encounter errors.  When you do, please indicate what's
% wrong and how to fix it

a = 1;

a

a = 1

a = [1:3:100]

whos a

a'

whos a

a = a/10

b = [1:3]

a/b

figure
plot(a)

figure
plot(sin(a))

figure
plot(a, sin(a))

b = 3 * ones(size(a))

c = a .* b;

whos c

figure
plot(c)

c = a * b;

whos a b

c = a' * b ;

whos c

c = b' * a ;

whos c

figure
imagesc(c)

string = 'SPM_fMRIDesMtx'

whos string

string = [string '.mat']

path_string = pwd

string = sprintf('%s/%s', path_string, string)

help sprintf

string = sprintf('vol_blahblahblah_%04d.img', 45)

who

clear

whos

load SPM_fMRIDesMtx

whos

X = xX.X


figure
imagesc(X)

r1 = X(:,1); 

r2 = X(:,2);

figure
plot(r1,'r');
hold on
plot(r2, 'b');
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 2 : creating a model for the BOLD signal corresponding to %
% a train of stimuli

figure
% Our resolution will be 100 points per second.
points = 100;

% create the HRF using SPM's function
help spm_hrf

% Fill in the blank  (the ???? , that is):
hrf = spm_hrf( 1/points );
subplot 311, plot(hrf),title('Response to a single stimulus')
axis([0 5000 -2e-3 4e-3])

% create a set of stimuli occurring at random times 
% in the time series
stim=zeros(50*points,1);

% Put ones where you want to have stimuli occur
stim( 5*points:13*points:end ) =1;

subplot 312, plot(stim),title('A set of stimuli at randomized times')


% convolve the input and the HRF and clip out the end (empty)
help conv
resp = conv(stim, hrf);

resp = resp(1:length(stim));
subplot 313, plot(resp),title('Response to the set of stimuli')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 3 : reading in a time series of images from a single file and 
% displaying the images

help read_hdr

hdr = read_hdr('tmp.hdr')

xdim = hdr.xdim;
ydim = hdr.ydim;
zdim = hdr.zdim;
tdim = hdr.tdim;

whos timeseries

pF = fopen('tmp.img','rb','ieee-be');
timeseries = (fread(pF,[xdim*ydim*zdim*tdim], 'short'))'; 	
fclose(pF);

whos timeseries

if tdim >=2
    timeseries = reshape(timeseries, xdim*ydim*zdim, tdim);
    timeseries = timeseries';
end

whos timeseries

help sub2ind
index = sub2ind([xdim ydim zdim], 28, 10, 25)

figure
plot(timeseries(:,index ))

title('time series at pixel 28,10,25')

oneFrame = timeseries( 5,:);
whos oneFrame

oneFrame = reshape(oneFrame,xdim,ydim,zdim );

whos oneFrame

figure
imagesc(squeeze(oneFrame(:,:,20)))
colormap(gray)

% Let's make a movie
for slice=1:zdim
    imagesc(squeeze(oneFrame(:,:,slice)))
    title(sprintf('slice %d', slice))
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3B : writing  a function
help function

orthogonal('tmp',32,32,14)

function result = orthogonal(root, x, y, z)
% function result = orthogonal(root, x, y, z)
%
% this function takes a file in analyze format and 
% displays orthogonal sections thru the planes interscting
% at x,y,z
% 
% root:  the name of the image file you want without the .img
% x, y, z : the coordinates of the point at which the three planes
% intersect.

	%figure out the whole name of the header and image files
	% and read them.  Note the use of strcat
	 h = read_hdr(strcat(root,'.hdr'));
     d = read_img2(h, (strcat(root,'.img') ));

    colormap(gray)
    
    % display the orthogonal sections
    subplot 221, imagesc(squeeze(d(:,:,z))')
    subplot 222, imagesc(squeeze(d(:,y,:))')
    subplot 223, imagesc(squeeze(d(x,:,:))')

    
return


