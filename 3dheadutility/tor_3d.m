function [D,Ds,hdr,p,coords,X,Y,Z] = tor_3d(varargin)
% function [D,Ds,hdr,p,coords,X,Y,Z] = tor_3d(varargin)
% made to use single_subj_T1.img from SPM99
%
% Made by Tor Wager, 10/3/2001 last modified 10/19/01
%
% options:
% 'data'        followed by image data, must also use 'hdr'; data is full image volume
% 'hdr'         followed by hdr structure (use read_hdr)
% 'coords'      3 element row vector of x,y,z coordinates at which to cut away
% 'figure'      create a new figure to plot on
% 'whichcuts'   followed by incisions; choices are x y z w (whole)
%               example:'whichcuts','xyz' (xyz is default)
%               order of output handles (p(i)) is 'wyzx'
% 'filename'    followed by filename of image data to load, in single quotes
%               should be analyze img format, without the .img extension
%               cluster imaging assumes neurological orientation, but should work anyway.
%
% output:
%   D = img data, Ds = smoothed data, hdr = header, p = image handles, coords = coordinates
%   X,Y,Z are the millimeter, origin centered reference frame for the head isosurface
%
% examples:
% [D,Ds,hdr,p,coords] = tor_3d('figure','data',D,'hdr',hdr,'whichcuts','yzx');
% [D,Ds,hdr,p,coords] = tor_3d('figure');

%------------------------------------------------------------------------------
% Set up arguments and default values
%------------------------------------------------------------------------------
D = [];
Ds =[];
coords = [0 0 0];
whichcuts = 'xyz';
filename = 'single_subj_T1';    % default structural image to use
topmm = 60;                     % image only this high in mm to avoid image distortions above head.

for i = 1:nargin
    if isstr(varargin{i})
        switch varargin{i}
        case 'data', D = varargin{i+1};
        case 'hdr', hdr = varargin{i+1};
        case 'coords', coords = varargin{i+1};
        case 'figure', h3dfig = figure; set(h3dfig,'Tag','myFig');,
            colormap(gray(100)),set(gcf,'Color','k'),axis off,set(gcf,'Position',[184   115   672   543])
        case 'whichcuts',whichcuts = varargin{i+1};
        case 'filename',filename = varargin{i+1};
        end
    end
end

%------------------------------------------------------------------------------
% Load the image file, if necessary
%------------------------------------------------------------------------------
if isempty(D)
    fullpath = which([filename '.img']);
    if isempty(fullpath), error(['Cannot find file: ' filename '.img']);
    else disp(['Loading structural image: ' fullpath]);drawnow
    end
    [array,hdr] = readim2(filename);
    % rotate 270 degrees so that front of head is positive y
    % (works with canonical SPM images, at least, so this orientation is ok.)
    for i = 1:size(array,3)
        D(:,:,i) = rot90(rot90(rot90((array(:,:,i)))));
    end

end

if isempty(Ds)
    Ds = smooth3(D);
end


%------------------------------------------------------------------------------
% Define X Y Z coordinates of head data
%------------------------------------------------------------------------------
% Define X, Y, Z relative to the origin, so head will be centered on origin
% Multiply by voxel size in header so that coordinates are in mm from the origin
% Make sure origin is set properly in header for this to work accurately.
[M N P] = size(D);
[X Y Z] = meshgrid(((1:N)-hdr.origin(1))*hdr.xsize, ((1:M)-hdr.origin(2))*hdr.ysize, ((1:P)-hdr.origin(3))*hdr.zsize);

%------------------------------------------------------------------------------
% Make patches for cutaway
%------------------------------------------------------------------------------
index = 1;
if any(whichcuts == 'w')
    [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[nan nan nan nan nan topmm]);          % whole head
    index = index + 2;
end

% xyz inset cutaway - seems to work only if x does not come first
if any(whichcuts == 'y')
    [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[nan nan nan coords(2) nan topmm]);
    index = index + 2;
end
 if any(whichcuts == 'z')
    [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[nan nan nan nan nan coords(3)]);
    index = index + 2;
end
if any(whichcuts == 'x')
    [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[nan coords(1) nan nan nan topmm]);
    index = index + 2;
end   


%------------------------------------------------------------------------------
% Set lighting conditions
%------------------------------------------------------------------------------

myLight = standardMRIlighting('full',p(1:2));
if length(p) > 2, standardMRIlighting('reflectance',p(3:4));, end
if length(p) > 4, standardMRIlighting('reflectance',p(5:6));, end

set(myLight,'Tag','myLight')

rotate3d on
	%------------------------------------------------------------------------------
    % set callback to light follow the camera
	%------------------------------------------------------------------------------
	if exist('lightFollowView') == 2
        set(gcf, 'WindowButtonUpFcn', 'lightFollowView');
    else
        warning('Cannot find lightFollowView.m to set light position.')
    end


%------------------------------------------------------------------------------
% Sub-function for imaging patch
%------------------------------------------------------------------------------
function [p1,p2] = imagePatch(X,Y,Z,D,Ds,inValues)

    headColor = [.8,.5,.4]; % old: [1 .75 .65]
    surface_threshold = mean(mean(mean(D))); 

    if isempty(inValues)
        inValues = [nan nan nan nan nan nan];
    end
    
    if isempty(X) | isempty(Y) | isempty(Z)
        % no xyz coordinates
        FV2 = isosurface(Ds,surface_threshold);
        IS2 = isocaps(D,surface_threshold);

        % Draw figure patches
	try
        	p1 = patch(FV2,'FaceColor',headColor,'EdgeColor','none','FaceAlpha',1,'SpecularExponent',200,'SpecularStrength',.2);
        	p2 = patch(IS2,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
	catch
        	p1 = patch(FV2,'FaceColor',headColor,'EdgeColor','none','SpecularExponent',200,'SpecularStrength',.2);
        	p2 = patch(IS2,'FaceColor','interp','EdgeColor','none');
	end
    else
        
        % Define subvolume for cutaway view
        [x y z A] = subvolume(X,Y,Z,D,inValues);
        [x y z As] = subvolume(X,Y,Z,Ds,inValues);
        FV2 = isosurface(x,y,z,As,surface_threshold);
        IS2 = isocaps(x,y,z,A,surface_threshold);

        % Draw figure patches
	try
        	p1 = patch(FV2,'FaceColor',headColor,'EdgeColor','none','FaceAlpha',1,'SpecularExponent',200,'SpecularStrength',.2);
        	p2 = patch(IS2,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
	catch
        	p1 = patch(FV2,'FaceColor',headColor,'EdgeColor','none','FaceAlpha',1,'SpecularExponent',200,'SpecularStrength',.2);
        	p2 = patch(IS2,'FaceColor','interp','EdgeColor','none');
	end
        % isonormals(A,p1)
    end
    
    drawnow
return



