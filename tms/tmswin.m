function fig = tmswin2()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.

load tmswin2

h0 = figure('Color',[0 0 0], ...
	'Colormap',mat0, ...
	'CreateFcn','tmscall(''create'');tmsglbls', ...
	'Position',[224 216 971 632], ...
	'Renderer','zbuffer', ...
	'RendererMode','manual', ...
	'Tag','Fig1');
h1 = uimenu('Parent',h0, ...
	'Label','TMS Files', ...
	'Tag','filetag');
h2 = uimenu('Parent',h1, ...
	'Callback','tmscall(''load_tms'')', ...
	'Label','Load Coordinates and Intensities', ...
	'Tag','LoadTMSTag');
h2 = uimenu('Parent',h1, ...
	'Callback','tmscall(''load_xyz'');', ...
	'Label','Load Polhemus Coordinates', ...
	'Tag','LoadXYZTag');
h2 = uimenu('Parent',h1, ...
	'Callback','tmscall(''save_file'');', ...
	'Label','Save Coordinates and Intensity', ...
	'Tag','saveTag');
h2 = uimenu('Parent',h1, ...
	'Callback','tmscall(''exit'');', ...
	'Label','Exit', ...
	'Tag','exitTag');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.5 0 0], ...
	'ListboxTop',0, ...
	'Position',[503.2 202.4 212 60.8], ...
	'Style','frame', ...
	'Tag','Frame1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','rotax(2);', ...
	'ListboxTop',0, ...
	'Max',360, ...
	'Position',[9.6 32 14.4 379.2], ...
	'Style','slider', ...
	'Tag','Slider2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','rotax(1)', ...
	'ListboxTop',0, ...
	'Max',360, ...
	'Position',[74.40000000000001 456.8 358.4 16.8], ...
	'Style','slider', ...
	'Tag','Slider1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'ButtonDownFcn','tmscall(''ClearFiles'');', ...
	'Callback','tmscall(''ClearFiles'')', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 69.59999999999999 132.8 25.6], ...
	'String','Clear', ...
	'Tag','Pushbutton1');
h1 = axes('Parent',h0, ...
	'Units','points', ...
	'View',[0 -1.8252], ...
	'CameraUpVector',[0 0 1], ...
	'Color',[0 0 0], ...
	'ColorOrder',mat1, ...
	'DataAspectRatioMode','manual', ...
	'NextPlot','add', ...
	'PlotBoxAspectRatioMode','manual', ...
	'Position',[73.59999999999999 43.2 375.2 386.4], ...
	'Tag','Axes1', ...
	'WarpToFill','off', ...
	'WarpToFillMode','manual', ...
	'XColor',[1 0 0], ...
	'YColor',[1 0 0], ...
	'ZColor',[1 0 0]);
h2 = text('Parent',h1, ...
	'Color',[1 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[0.4978602842106024 -8.137186947187368 -0.3618169010579003], ...
	'String','X', ...
	'Tag','Text1', ...
	'VerticalAlignment','cap');
set(get(h2,'Parent'),'XLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[1 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',mat2, ...
	'Rotation',-90, ...
	'String','Y', ...
	'Tag','Text2', ...
	'VerticalAlignment','baseline');
set(get(h2,'Parent'),'YLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[1 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[-0.08200269471614738 -8.155315039347878 0.2070587325202504], ...
	'Rotation',90, ...
	'String','Z', ...
	'Tag','Text3', ...
	'VerticalAlignment','baseline');
set(get(h2,'Parent'),'ZLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[0.4978602842106024 -8.172693473561898 0.7524094338827632], ...
	'Tag','Text4', ...
	'VerticalAlignment','bottom');
set(get(h2,'Parent'),'Title',h2);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','tmscall(''remove_junk'');', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 418.4 101.6 19.2], ...
	'String','Remove Points', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','tmscall(''correct_motion'');', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 442.4 100.8 20.8], ...
	'String','Motion Correction', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','tmscall(''fit'');', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 367.2 100.8 20.8], ...
	'String','Corregistration', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','tmscall(''get_area'');', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[505.6 238.4 100.8 21.6], ...
	'String','Surface Area', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','tmscall(''get_center'');', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 212.8 100.8 20.8], ...
	'String',' Activation Center', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','tmscall(''get_intensities'');', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 393.6 100.8 20.8], ...
	'String','Add Intensity File', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'ButtonDownFcn','tmscall(''ChangeActiveFile'');get(gcbo,''Value'')', ...
	'Enable','inactive', ...
	'Interruptible','off', ...
	'Max',5, ...
	'Min',1, ...
	'Position',[506.4 104.8 131.2 62.4], ...
	'Style','listbox', ...
	'Tag','Listbox1', ...
	'Value',2);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0 0 0], ...
	'ForegroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 170.4 127.2 20], ...
	'String','Active File', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','tmscall(''project'')', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[506.4 340.8 100.8 20.8], ...
	'String','Projection', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'ListboxTop',0, ...
	'Position',[506.4 264 100.8 21.6], ...
	'String','not used', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.7019609999999999 0.7019609999999999 0.7019609999999999], ...
	'Callback','tmscall(''make_analyze'')', ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[505.6 314.4 100.8 20.8], ...
	'String','Convert to Analyze', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Enable','inactive', ...
	'ListboxTop',0, ...
	'Position',[506.4 290.4 100.8 20], ...
	'String','not used', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0 0 0.1], ...
	'FontSize',12, ...
	'ForegroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[617.6 208.8 91.2 49.6], ...
	'Style','text', ...
	'Tag','ResultBox');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'ForegroundColor',[0 0 1], ...
	'ListboxTop',0, ...
	'Position',[612.8000000000001 342.4 88.80000000000001 16], ...
	'String','Expansion', ...
	'Style','checkbox', ...
	'Tag','ExpansionCheckBox');
if nargout > 0, fig = h0; end
