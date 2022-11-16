% This ui is set to render on the default single_subj_T1, in default colors
% More flexibility is available if you use the functions.
% By Tor Wager, 10/3/2001, last edit 10/22/01
%
% main functions used:
%   tor_3d.m - images head with cutaway views
%   imageCluster.m - images a cluster isosurface
%   mni_TSU.m and tor_ihb_TalSpace.m - to get clusters
%       these and related functions are part of Talairach Space Utility
%       written by Sergey Pakhomov, 2001
%       modified very slightly by Tor Wager to not convert to Talairach Space
%       and use MNI coordinates instead.
%
%   Use of the TSU functions for getting clusters require SPM99.
%
% Output to workspace:
%   Isosurface handles are in p, for head isosurfaces, and cH, for the cluster isosurface
%   D = image data, Ds = smoothed data, hdr = img header
%
% Add more clusters by using:
% cH(2) = imageCluster('cluster',clusters(i));

clear coords
clear cl_center
myColors = {'r' 'r' 'r' 'r' 'b' 'b' 'b'};


%------------------------------------------------------------------------------
% Get head, if not default
%------------------------------------------------------------------------------
cflag = 0;
while cflag == 0
    a = input('Use default head single_subj_T1...?  (type default or file name, no .img extension)','s');
    switch a
    case 'default', fileName = 'single_subj_T1';
    otherwise fileName = a;
    end
    if exist([fileName '.img']) == 2, 
	cflag = 1;
    else
	disp(['Can''t locate file: ' fileName '.img'])
    end
end


%------------------------------------------------------------------------------
% Add surface?
%------------------------------------------------------------------------------
cflag = 0;
while cflag == 0
    a = input('Add brain surface (y\n\filename,no .img extension)...?','s');
    switch a
    case 'y', fileName2 = ['surf_' fileName];, dosurf = 1;
    case 'n', fileName2 = [];, cflag = 1;, dosurf = 0;
    otherwise fileName2 = a;, dosurf = 1;
    end
    if dosurf & exist([fileName2 '.mat']) == 2, 
	    cflag = 1;
    elseif dosurf
	    disp(['Can''t locate file: ' fileName '.mat'])
    end
end



%------------------------------------------------------------------------------
% Get clusters, if necessary
%------------------------------------------------------------------------------
cflag = 0;
while cflag == 0
    a = input('Get clusters from...?  (workspace,file,figure)','s');
    switch a
    case 'workspace', cflag = 1;
    case 'file', clusters = imageCluster('getclusters');, cflag = 1;
    case 'figure', clusters = imageCluster('getfigclusters');, cflag = 1;
    otherwise disp(['Not a valid choice...type entry string with no quotes.'])
    end
end

%------------------------------------------------------------------------------
% Get number of cluster to image 
%------------------------------------------------------------------------------
cluster_names = str2mat(clusters.name);
clNums = num2str((1:length(clusters))');
blankSpcs = repmat(' ',length(clusters),1);
for i = 1:length(clusters),cl_center{1,i} = num2str(round(mean(clusters(i).XYZmm')));,end
cl_center = str2mat(cl_center);
cluster_names_wNums = [clNums blankSpcs cluster_names blankSpcs cl_center]
clear clNums,clear blankSpcs

cflag =0;
while cflag == 0
    a = input('Enter number of cluster to image, or vector of clusters, ex. [2 4] : ');
    if any(a > length(clusters)), disp(['Invalid: no cluster with number ' num2str(a(a>length(clusters))) '.']),
    else cflag = 1;     % ok to go ahead
    end
end


%------------------------------------------------------------------------------
% Define cluster and cluster center (-5mm for easy viewing)
%------------------------------------------------------------------------------
clOrig = clusters(a);
cl = [];

disp(['Selected ' num2str(length(a)) ' clusters.'])
for i = 1:length(a)

    %------------------------------------------------------------------------------
    % get xMin xMax and vTal if not already defined (if clusters made with RoiUtility)
    %------------------------------------------------------------------------------
    if ~(isfield(clOrig(i),'xMin'))
	    mycl = get_cluster_volume(clOrig(i));
    else
	    mycl = clOrig(i);
    end
    if i == 1,
	    cl = mycl;
    else
	    cl(i) = mycl;
    end 

    coords(i,1) = mean([cl(i).xMin cl(i).xMax]) - 5;
    coords(i,2) = mean([cl(i).yMin cl(i).yMax]) - 5;
    coords(i,3) = mean([cl(i).zMin cl(i).zMax]) - 5;
end

disp(['Cluster centers are:'])
coords
if size(coords,1) > 1
    bestCoords = min(coords);
else
    bestCoords = coords;
end
disp(['Coordinates to use are: ' num2str(bestCoords)])

%------------------------------------------------------------------------------
% Get views to cut away
%------------------------------------------------------------------------------
whichc = input('Enter axes to cut along.  Valid choices are xyzw. Example yz: ','s');

%------------------------------------------------------------------------------
% Image the head and add the cluster
%------------------------------------------------------------------------------
% set clusters first or make3davi won't work. ??? or maybe it's the figure size.
cH = [];
for i = 1:length(cl)
    cH(i) = imageCluster('cluster',cl(i),'color',myColors{i},'alpha',.7); drawnow
    set(cH(i),'Tag',['cluster' num2str(i)])
end
[D,Ds,hdr,p,bestCoords] = tor_3d('whichcuts',whichc,'coords',bestCoords,'filename',fileName);
colormap(gray(100))

O.H = cH;

if dosurf
    eval(['load ' fileName2])
    try
        p(end+1) = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5],'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',.6,'SpecularExponent',200);
    catch
        p(end+1) = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5],'EdgeColor','none','SpecularStrength',.2,'SpecularExponent',200);
    end
end  

saveas(gcf,'temp','fig')
save temp_handles p cH

%set(p(1),'SpecularExponent',200)
%set(p(3),'SpecularExponent',200)
%set(p(5),'SpecularExponent',200)
axis off
set(gcf,'Color','k')
%set(p(1),'FaceColor',[.8 .5 .4]);
%set(p(3),'FaceColor',[.8 .5 .4]);
%set(p(5),'FaceColor',[.8 .5 .4]);
%set(p(1),'SpecularStrength',.2)
%set(p(3),'SpecularStrength',.2)
%set(p(5),'SpecularStrength',.2)

% make3davi(O)