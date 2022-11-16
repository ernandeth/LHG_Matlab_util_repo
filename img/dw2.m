function fig = dw2()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.

load dw2

h0 = figure('Color',[0 0 0], ...
	'Colormap',mat0, ...
	'FileName','/net/pabst/y/hernan/matlab/img/dw2.m', ...
	'PaperPosition',[-22.02352941176471 88.94117647058823 656.4705882352941 614.964705882353], ...
	'PaperUnits','points', ...
	'Position',[352 132 775 702], ...
	'Tag','MainDispWinTag', ...
	'ToolBar','none');
h1 = uimenu('Parent',h0, ...
	'Label','DispWinFile', ...
	'Tag','filetag');
h2 = uimenu('Parent',h1, ...
	'Callback','dwcall(''load_*img'')', ...
	'Label','Load *img', ...
	'Tag','DisplayTag');
h2 = uimenu('Parent',h1, ...
	'Label','Load SPM*.mat', ...
	'Tag','CreateSpmTag');
h3 = uimenu('Parent',h2, ...
	'Callback','dwcall(''load_spmF'')', ...
	'Label','F statistic', ...
	'Tag','loadF');
h3 = uimenu('Parent',h2, ...
	'Callback','dwcall(''load_spmt'')', ...
	'Label','t statistic', ...
	'Tag','loadT');
h2 = uimenu('Parent',h1, ...
	'Callback','dwcall(''save'')', ...
	'Label','Save *img', ...
	'Tag','saveTag');
h2 = uimenu('Parent',h1, ...
	'Callback','dwcall(''exit'');', ...
	'Label','Exit', ...
	'Tag','ExitTag');
h1 = uimenu('Parent',h0, ...
	'Label','ROI', ...
	'Tag','roitag');
h2 = uimenu('Parent',h1, ...
	'Label','Polygon', ...
	'Tag','ROIPolygon1');
h2 = uimenu('Parent',h1, ...
	'Label','Oval', ...
	'Tag','ROIOval1');
h2 = uimenu('Parent',h1, ...
	'Label','Square', ...
	'Tag','ROISquare1');
h1 = uimenu('Parent',h0, ...
	'Label','Flip', ...
	'Tag','FlipTag');
h2 = uimenu('Parent',h1, ...
	'Label','Flip Z', ...
	'Tag','FlipZTag');
h2 = uimenu('Parent',h1, ...
	'Label','Flip Y', ...
	'Tag','flipYTag');
h2 = uimenu('Parent',h1, ...
	'Label','Flip X', ...
	'Tag','FlipXTag');
h2 = uimenu('Parent',h1, ...
	'Label','uimenu', ...
	'Tag','Flipuimenu1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'ListboxTop',0, ...
	'Position',[31.2 77.59999999999999 12.8 276.8], ...
	'Style','slider', ...
	'Tag','Slider1', ...
	'Value',0.002);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'ListboxTop',0, ...
	'Position',[14.4 77.59999999999999 12.8 276.8], ...
	'Style','slider', ...
	'Tag','Slider1', ...
	'Value',0.1);
h1 = axes('Parent',h0, ...
	'Units','points', ...
	'Box','on', ...
	'CameraUpVector',[0 -1 0], ...
	'CameraUpVectorMode','manual', ...
	'Color',[1 1 1], ...
	'ColorOrder',mat1, ...
	'Layer','top', ...
	'NextPlot','add', ...
	'Position',[64 73.59999999999999 416 424.8], ...
	'Tag','Axes1', ...
	'XColor',[0 0 0], ...
	'XLim',[0.5 512.5], ...
	'XLimMode','manual', ...
	'YColor',[0 0 0], ...
	'YDir','reverse', ...
	'YLim',[0.5 512.5], ...
	'YLimMode','manual', ...
	'ZColor',[0 0 0]);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[255.4551020408164 532.98 9.160254037844386], ...
	'Tag','Axes1Text4', ...
	'VerticalAlignment','cap');
set(get(h2,'Parent'),'XLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[-27.71224489795919 257.5239999999999 9.160254037844386], ...
	'Rotation',90, ...
	'Tag','Axes1Text3', ...
	'VerticalAlignment','baseline');
set(get(h2,'Parent'),'YLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','right', ...
	'Position',[-79.95714285714287 -115.2120000000001 9.160254037844386], ...
	'Tag','Axes1Text2', ...
	'Visible','off');
set(get(h2,'Parent'),'ZLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[255.4551020408164 -4.620000000000005 9.160254037844386], ...
	'Tag','Axes1Text1', ...
	'VerticalAlignment','bottom');
set(get(h2,'Parent'),'Title',h2);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'Callback','dwcall(''addspm'')', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 446.4 91.2 20.8], ...
	'String','Add SPM', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','dwcall(''pixel_count'')', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 411.2 91.2 20.8], ...
	'String','Count Pixels ', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','dwcall(''multislice_toggle'')', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[505.6 378.4 91.2 20.8], ...
	'String','Multi Slice', ...
	'Style','checkbox', ...
	'Tag','msToggle', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'Callback','dwcall(''get_coords'')', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 345.6 91.2 20.8], ...
	'String','Get Coordinates', ...
	'Style','checkbox', ...
	'Tag','CoordsButton');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','dwcall(''reset'')', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 310.4 91.2 20.8], ...
	'String','Reset', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'Callback','dwcall(''timeseries'');', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[503.2 478.4 92 20], ...
	'String','Time Plot', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','dwcall(''changeslice'')', ...
	'ForegroundColor',[1 0 0], ...
	'ListboxTop',0, ...
	'Max',60, ...
	'Min',1, ...
	'Position',[520.9411764705882 221.9294117647059 55.90588235294118 20.32941176470588], ...
	'String',['3';' '], ...
	'Style','edit', ...
	'Tag','EditSliceNum');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','dwcall(''back'')', ...
	'ListboxTop',0, ...
	'Position',[489.6 220 31.2 22.4], ...
	'String','<<', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','dwcall(''forward'')', ...
	'ListboxTop',0, ...
	'Position',[576.8 220.8 31.2 22.4], ...
	'String','>>', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0 0 0], ...
	'ForegroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[45.6 520.8 443.2 15.2], ...
	'Style','text', ...
	'Tag','ImageFileName');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0 0 0], ...
	'ForegroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[50.4 35.2 443.2 15.2], ...
	'Style','text', ...
	'Tag','SPMFileName');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ForegroundColor',[1 0 0], ...
	'ListboxTop',0, ...
	'Max',60, ...
	'Min',1, ...
	'Position',mat2, ...
	'String','32.8259      50.5758', ...
	'Style','text', ...
	'Tag','pixel_coordinates');
if nargout > 0, fig = h0; end