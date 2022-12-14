function tms3d()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.

load disp3d                             

a = figure('Color',[1 1 1], ...
	'Colormap',mat0, ...
	'Position',[584 294 835 555], ...
	'Renderer','zbuffer', ...
	'RendererMode','manual', ...
	'Tag','Fig1');
b = uicontrol('Parent',a, ...
	'Units','points', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'Callback','rotax(2);', ...
	'Max',180, ...
	'Min',-180, ...
	'Position',[469.6 359.2 150.4 19.2], ...
	'Style','slider', ...
	'Tag','Slider2', ...
	'Value',-50.454);
b = uicontrol('Parent',a, ...
	'Units','points', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'Callback','rotax(1)', ...
	'Max',180, ...
	'Min',-180, ...
	'Position',[468.8 398.4 150.4 20], ...
	'Style','slider', ...
	'Tag','Slider1', ...
	'Value',-100.829);
b = uicontrol('Parent',a, ...
	'Units','points', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'Callback','clear all;close;', ...
	'Position',[468 44.8 154.4 40.8], ...
	'String','CLOSE', ...
	'Tag','Pushbutton1');
b = axes('Parent',a, ...
	'Units','points', ...
	'CameraUpVector',[0 0 1], ...
	'CLimMode','manual', ...
	'Color',[1 1 1], ...
	'ColorOrder',mat1, ...
	'Position',[47.2 36.8 398.4 393.6], ...
	'Tag','Axes1', ...
	'View',[322.5 30], ...
	'XColor',[0 0 0], ...
	'XLim',[-15 10], ...
	'XLimMode','manual', ...
	'YColor',[0 0 0], ...
	'YLim',[-10 10], ...
	'YLimMode','manual', ...
	'ZColor',[0 0 0], ...
	'ZLim',[0 20], ...
	'ZLimMode','manual');
c = line('Parent',b, ...
	'Color',[0 0 1], ...
	'LineStyle','none', ...
	'Marker','o', ...
	'Tag','Axes1Line1', ...
	'XData',mat2, ...
	'YData',mat3, ...
	'ZData',mat4);
c = text('Parent',b, ...
	'Color',[0 0 0], ...
	'HandleVisibility','callback', ...
	'Position',[-8.86487 -15.8581 0], ...
	'String','X', ...
	'Tag','Axes1Text4', ...
	'VerticalAlignment','top');
set(get(c,'Parent'),'XLabel',c);
c = text('Parent',b, ...
	'Color',[0 0 0], ...
	'HandleVisibility','callback', ...
	'HorizontalAlignment','right', ...
	'Position',[-22.2294 -6.52767 0], ...
	'String','Y', ...
	'Tag','Axes1Text3', ...
	'VerticalAlignment','top');
set(get(c,'Parent'),'YLabel',c);
c = text('Parent',b, ...
	'Color',[0 0 0], ...
	'HandleVisibility','callback', ...
	'HorizontalAlignment','center', ...
	'Position',[-3.59174 24.6746 0], ...
	'Rotation',90, ...
	'String','Z', ...
	'Tag','Axes1Text2', ...
	'VerticalAlignment','baseline');
set(get(c,'Parent'),'ZLabel',c);
c = text('Parent',b, ...
	'Color',[0 0 0], ...
	'HandleVisibility','callback', ...
	'HorizontalAlignment','center', ...
	'Position',[35.0845 39.2312 0], ...
	'Tag','Axes1Text1', ...
	'VerticalAlignment','bottom');
set(get(c,'Parent'),'Title',c);
