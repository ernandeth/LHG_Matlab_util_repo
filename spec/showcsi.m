% showcsi(data)
% - For visualizing 4D data sets.
% - Allows basic manipulations such as FFT, cropping,
%   zero-filling and cosine windowing.
%
% (C) 2000 Jeffrey Tsao
% Biomedical Magnetic Resonance Laboratory
% University of Illinois at Urbana-Champaign
% jtsao2@hotmail.com

% MODIFICATION HISTORY
% Aug 31, 2000
%  - Fixed bug in handling Z dimension
%  - Fixed bug in handling size of 1 for X & Y dimensions
% Jan 26, 2001
%  - Fixed bug in cropping -'ve
%  - Made the code more or less self-sufficient by including
%    the cosine, transformKspaceToImage, transformImageToKspace,
%    shift functions
% Feb 8, 2001
%  - Fixed bug: didn't display properly when switching modes
%    for complex values e.g. abs -> phase
%    Needed to set image axes as current axes in SetUpColorScale
% Aug 13, 2001
%  - Removed all findobj tag statements so that multiple instances of
%    this program can run at the same time
%  - Added the "Close All" button
% Feb 13, 2002
%  - Fixed bug so that when you change frequency position, the
%    slice is to change back to center slice. SetUpImage uses
%    the current slice (ud.PosZ) instead.
%  - Change frequency offset to refer to center frequency instead
%  - Change frequency offset so that it is according to sign of the spectral width
%  - Change the export variable name to match input variable name, if there is one.
%  - Added 'Parent' attribute when calling uicontrol so that interface elements are
%    always created in the correct window.
% Feb 19, 2002
%  - Fixed zeropadding in SetUpImage. Had error if one of the padding size had a size of 0 before.
% Feb 22, 2002
%  - Switch Hamming to sqrt(Hamming) filtering
% Feb 28, 2002
%  - Update spectral bandwidth and center position when cropped.
%    Cropping direction fixed depending on sign of bandwidth.s

function showcsi(action, varargin)
  if nargin<1,
    help(mfilename)
    return
  end
  if isequal(class(action), 'char'),  % A command
    feval(action,varargin{:})
  else
    if ndims(action)>4,  % action is a variable
      error('Cannot handle data with more than 4 dimensions (3 spatial + 1 spectral)');
    end
    Initialize(action, inputname(1));
  end
  clear action
return;
%####################################################################
%####################################################################
%####################################################################
function Initialize(data, dataname)
  if isempty(dataname), titleString='CSI'; else titleString=['CSI - ' dataname]; end 
  
  % If program is already running, bring it to the foreground
  %h = findobj(allchild(0), 'tag', 'CSI');
  %if isempty(h)
  %  MainFig = figure('tag', 'CSI','IntegerHandle', 'off'); % New figure
  %else
  %  MainFig = figure(h(1));  % Pop it up
  %  delete(get(MainFig,'children'));
  %end
  %clear h
  
  MainFig = figure('tag', 'CSIwin','IntegerHandle', 'off'); % Always open a new figure
  if get(0, 'ScreenDepth')>8  grayres=256;
                         else grayres=128; end
  set(MainFig, ...
      'Name',titleString, ...
      'NumberTitle','off', 'HandleVisibility', 'on', ...
      'Visible','off', 'Resize', 'off',...
      'BusyAction','Queue','Interruptible','off', ...
      'Color', [.8 .8 .8], ...
      'MenuBar','figure', ...  % disable tools menu
      'IntegerHandle', 'off', ...
      'Colormap', gray(grayres),...
      'WindowButtonMotionFcn', [mfilename '( ''TrackMouseMotion'')'],...
      'WindowButtonDownFcn', [mfilename '( ''TrackMouseButton'')'] );         
  clear grayres titleString
    
  %---- SET UP FIGURE ------------------------------------
  FigColor = get(MainFig, 'Color');    
  figpos = get(MainFig, 'position');
  figpos(3:4) = [560 420];
  % Adjust the size of the figure window
  horizDecorations = 10;  % resize controls, etc.
  vertDecorations = 45;   % title bar, etc.
  screenSize = get(0,'ScreenSize');
  if (screenSize(3) <= 1)
    % No display connected (apparently)
    screenSize(3:4) = [100000 100000]; % don't use Inf because of vms
  end
  if (((figpos(3) + horizDecorations) > screenSize(3)) | ...
      ((figpos(4) + vertDecorations) > screenSize(4)))
    % Screen size is too small for this window!
    delete(fig);
    error(['Screen resolution is too low ', ...
                '(or text fonts are too big) to run this program']);
  end
  dx = screenSize(3) - figpos(1) - figpos(3) - horizDecorations;
  dy = screenSize(4) - figpos(2) - figpos(4) - vertDecorations;
  if (dx < 0) figpos(1) = max(5,figpos(1) + dx); end
  if (dy < 0) figpos(2) = max(5,figpos(2) + dy); end
  set(MainFig, 'position', figpos);
  %  background
  clear screenSize figpos dx dy horizDecorations vertDecorations
  %------------------------------------------------------------------
  
  %------------------------------------------------------------------
  % Parameters for all buttons and menus
  Std.Interruptible = 'off';
  Std.BusyAction = 'queue';
  Std.Parent = MainFig;
  Btn = Std;
  Btn.Units = 'Pixels';
  Btn.Parent = MainFig;
  Btn.Style = 'pushbutton';
  Btn.Enable = 'on';
  Edt = Std;
  Edt.Units = 'Pixels';
  Edt.Parent = MainFig;
  Edt.Style = 'Edit';
  Edt.Enable = 'on';
  Rad = Btn;
  Rad.Style = 'radiobutton';
  %------------------------------------------------------------------
  
  %------------------------------------------------------------------
  % PROCESSING CONTROLS
  uicontrol(Btn, ...
    'Position',[2 365 35 20], ...
    'String','FFT', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''ProcessData'', ''FFT'')'] );
  uicontrol(Btn, ...
    'Position',[2 345 35 20], ...
    'String','IFFT', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''ProcessData'', ''IFFT'')'] );
  uicontrol(Btn, ...
    'Position',[2+35 365 85 20], ...
    'String','Shift to center', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''ProcessData'', ''Shift to center'')'] );
  uicontrol(Btn, ...
    'Position',[2+35 345 85 20], ...
    'String','0-pad to center', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''ProcessData'', ''Zeropad to center'')'] );  
  uicontrol(Btn, ...
    'Position',[2+120 365 80 20], ...
    'String','x Hamming^.5', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''ProcessData'', ''Multiply with Hamming'')'] );  
  uicontrol(Btn, ...
    'Position',[2+120 345 80 20], ...
    'String','/ Hamming^.5', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''ProcessData'', ''Divide by Hamming'')'] );  
  uicontrol(Btn, ...
    'Position',[2+200 365 60 20], ...
    'String','Crop +''ve', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''ProcessData'', ''Crop Positive'')'] );  
  uicontrol(Btn, ...
    'Position',[2+200 345 60 20], ...
    'String','Crop -''ve', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''ProcessData'', ''Crop Negative'')'] );  
  
  ud.ApplyToDimXh = uicontrol(Btn, ...
    'Style','Togglebutton', ...
    'Position',[2 385 20 20], ...
    'String','x', ...
    'Value',0, ...
    'Tag', '' );
  ud.ApplyToDimYh = uicontrol(Btn, ...
    'Style','Togglebutton', ...
    'Position',[2+20 385 20 20], ...
    'String','y', ...
    'Value',0, ...
    'Tag', '' );
  ud.ApplyToDimZh = uicontrol(Btn, ...
    'Style','Togglebutton', ...
    'Position',[2+40 385 20 20], ...
    'String','z', ...
    'Value',0, ...
    'Tag', '' );
  ud.ApplyToDimFh = uicontrol(Btn, ...
    'Style','Togglebutton', ...
    'Position',[2+60 385 20 20], ...
    'String','f', ...
    'Value',1, ...
    'Tag', '' );
  %------------------------------------------------------------------
  
  %------------------------------------------------------------------
  % IMAGE AND CONTROLS
  ud.ImageSliceCtrls_h = [];
  ud.ImageSliceCtrls_h(end+1) = uicontrol(Btn, ...
    'Position',[2 66 20 20], ...
    'String','|<', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''GoToSlice'', ''first'')'] );
  ud.ImageSliceCtrls_h(end+1) = uicontrol(Btn, ...
    'Position',[22 66 20 20], ...
    'String','<', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''GoToSlice'', ''prev'')'] );
  ud.ImageSliceCtrls_h(end+1) = uicontrol(Btn, ...
    'Position',[42 66 20 20], ...
    'String','>', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''GoToSlice'', ''next'')'] );
  ud.ImageSliceCtrls_h(end+1) = uicontrol(Btn, ...
    'Position',[62 66 20 20], ...
    'String','>|', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''GoToSlice'', ''last'')'] );
  ud.SliceNum_h = uicontrol(Std, ...
    'Parent', MainFig, ...
    'Style','text', ...
    'Units','pixels', ...
    'Position',[2 51 100 15], ...
    'Horiz','left', ...
    'Background', FigColor, ...
    'Tag', '', ...
    'String',{'' ''});
  ud.ImageControls_h = [];
  ud.ImageControls_h(end+1) = uicontrol(...
    'Parent', MainFig, ...
    'style','popup', ...
    'position',[156 66 54 20], ...
    'string',{'Gray';'Jet';'Hot';'Cool';'Bone';'Hsv'}, ...
    'Tag', '', ...
    'value',1, ...
    'callback',['tmpvar=get(gco,''string''); ' mfilename '( ''changeColors'', tmpvar{get(gco,''value'')} ); clear tmpvar']);
  ud.ImageComplex_h = uicontrol(...
    'Parent', MainFig, ...
    'style','popup', ...
    'position',[100 66 54 20], ...
    'string',{'Abs';'Real';'Imag.';'Phase'}, ...
    'Tag', '', ...
    'value',1, ...
    'Callback',[mfilename '( ''SetUpColorScale'',1 )'] );
  ud.ImageControls_h(end+1) = ud.ImageComplex_h;
  ud.ImageControls_h(end+1) = uicontrol(Btn, ...
    'Position',[210 66 26 20], ...
    'String','-''ve', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''invertColors'' )']);
  ud.ImageControls_h(end+1) = uicontrol(Btn, ...
    'Position',[238 76 20 10], ...
    'String','/\', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''changeBrightness'', 0.1 )']);
  ud.ImageControls_h(end+1) = uicontrol(Btn, ...
    'Position',[238 66 20 10], ...
    'String','\/', ...
    'Tag', '', ...
    'Callback',[mfilename '( ''changeBrightness'', -0.1 )']);
  %------------------------------------------------------------------

  %------------------------------------------------------------------
  % OVERALL INTERFACE ELEMENTS
  uicontrol(Btn, ...              % Help button
    'Position',[10 7 50 26], ...
    'String','Help', ...
    'Tag', 'ButHelp', ...
    'Callback',['helpwin ' mfilename]);
  uicontrol(Btn, ...              % Print button
    'Position',[60 7 55 26], ...
    'String','Print...', ...
    'Tag', 'ButHelp', ...
    'Callback','printdlg(''-crossplatform'',gcbf)')
%  uicontrol(Btn, ...              % Go to First Page button
%    'Position',[250 7 70 26], ...
%    'String','Start Over', ...
%    'Tag', 'ButFirstPage', ...
%    'Callback',[mfilename '( ''FirstPage'' )']);
%  uicontrol(Btn, ...              % Previous Page button
%    'Position',[330 7 70 26], ...
%    'String','< Previous', ...
%    'Tag', 'ButPrev', ...
%    'Callback',[mfilename '( ''PrevPage'' )']);
%  uicontrol(Btn, ...              % Next Page button
%    'Position',[400 7 70 26], ...
%    'String','Next >', ...
%    'Tag', 'ButNext', ...
%    'Callback',[mfilename '( ''NextPage'' )']);
  uicontrol(Btn, ...              % Exit button
    'Position',[500 7 50 26], ...
    'String','Close', ...
    'Tag', 'ButExit', ...
    'Callback','close(gcbf)');
   uicontrol(Btn, ...              % Exit button
    'Position',[430 7 70 26], ...
    'String','Close all...', ...
    'Tag', 'ButExit', ...
    'Callback','idx = findobj(allchild(0), ''tag'', ''CSIwin''); if length(idx)==1, close(idx); return; end; if isequal(questdlg(''Are you sure to close all plot windows?'',''Close all...'',''No''),''Yes''), close(findobj(allchild(0), ''tag'', ''CSIwin'')), end');
%  CSIwin
  %------------------------------------------------------------------
  uicontrol( Std, ...
      'Style','Text',...
      'Position',[2+300-5 345-256-45 90 20], ...
		'Horiz','right', ...
		'Background',FigColor, ...
        'Foreground','black', ...
        'String','Bandwidth (ppm)', ...
  		  'Horiz','left', ...
        'callback','');
  ud.InputSpecWidth_h = uicontrol( Edt, ...
      'Position',[2+300+90-5 345-256-45 50 20], ...
		'Horiz','right', ...
		'Background','white', ...
      'Foreground','black', ...
      'Tag','',...
      'callback', ['tmpvar=get(gco,''string''); ' mfilename '(''changeSpecWidth'', tmpvar); clear tmpvar']);
  uicontrol( Std, ...
      'Style','Text',...
      'Position',[2+300+90+50+15 345-256-45 45 20], ...
		'Horiz','right', ...
		'Background',FigColor, ...
        'Foreground','black', ...
        'String','Center', ...
  		  'Horiz','left', ...
        'callback','');
  ud.InputFreqOffset_h = uicontrol( Edt, ...
      'Position',[2+300+90+50+15+45 345-256-45 50 20], ...
		'Horiz','right', ...
		'Background','White', ...
      'Foreground','black', ...
      'Tag','',...
      'callback', ['tmpvar=get(gco,''string''); ' mfilename '(''changeFreqOffset'', tmpvar); clear tmpvar']);
  uicontrol(Btn, ...
    'Style','Togglebutton', ...
    'Position',[78+300-40 365+10 40 20], ...
    'String','Zoom', ...
    'ForegroundColor',[0,0,0],...
    'Value',0, ...
    'Tag', '',...
    'Callback', 'if get(gco,''Value''), zoom on; else zoom off; end');
    zoom off
  uicontrol(Btn, ...
    'Position',[78+300-90 365+10 50 20], ...
    'String','Save...', ...
    'ForegroundColor',[0,0,0],...
    'Callback', [mfilename '(''AskForVariableName'')'] );
   uicontrol(Btn, ...
    'Position',[78+300 365+10 30 20], ...
    'String','|<', ...
    'Tag', '', ...
    'Callback',     [mfilename '( ''ChangeFreqPos'', Inf)'] );
   uicontrol(Btn, ...
    'Position',[108+300 365+10 30 20], ...
    'String','<<', ...
    'Tag', '', ...
    'Callback',     [mfilename '( ''ChangeFreqPos'', 10)'] );
   uicontrol(Btn, ...
    'Position',[138+300 365+10 30 20], ...
    'String','<', ...
    'Tag', '', ...
    'Callback',     [mfilename '( ''ChangeFreqPos'', 1)'] );
   uicontrol(Btn, ...
    'Position',[168+300 365+10 30 20], ...
    'String','>', ...
    'Tag', '', ...
    'Callback',     [mfilename '( ''ChangeFreqPos'', -1)'] );
   uicontrol(Btn, ...
     'Position',[198+300 365+10 30 20], ...
    'String','>>', ...
    'Tag', '', ...
    'Callback',     [mfilename '( ''ChangeFreqPos'', -10)'] );
   uicontrol(Btn, ...
     'Position',[228+300 365+10 30 20], ...
    'String','>|', ...
    'Tag', '', ...
    'Callback',     [mfilename '( ''ChangeFreqPos'', -Inf)'] );
%------------------------------------------------------------------
  ud.PlotAbs_h = uicontrol(Btn, ...
    'Style','Togglebutton', ...
    'Position',[8+310 385+10 40 20], ...
    'String','Abs', ...
    'ForegroundColor',[0,0,0],...
    'Value',1, ...
    'Tag', '',...
    'Callback',     [mfilename '( ''SetUpGraph'' )'] );
  ud.PlotReal_h = uicontrol(Btn, ...
    'Style','Togglebutton', ...
    'Position',[8+350 385+10 40 20], ...
    'String','Real', ...
    'ForegroundColor',[0,0,1],...
    'Value',1, ...
    'Tag', '',...
    'Callback',     [mfilename '( ''SetUpGraph'' )'] );
  ud.PlotImag_h = uicontrol(Btn, ...
    'Style','Togglebutton', ...
    'Position',[8+390 385+10 40 20], ...
    'String','Imag', ...
    'ForegroundColor',[1,0,0],...    
    'Value',1, ...
    'Tag', '',...
    'Callback',     [mfilename '( ''SetUpGraph'' )'] );
%    'ForegroundColor',[0,0,0],...
%    'Color',[0,0,1],...
%    'Color',[1,0,0],...
  
  
  uicontrol(Btn, ...
    'Position',[258+180 385+10 60 20], ...
    'String','Scale X', ...
    'Callback', [mfilename '( ''AutoScaleAxis'', ''x'')']);
  uicontrol(Btn, ...
    'Position',[258+240 385+10 60 20], ...
    'String','Scale Y', ...
    'Callback', [mfilename '( ''AutoScaleAxis'', ''y'')']);
  
  %------------------------------------------------------------------
  
  %------------------------------------------------------------------
  % Set up filenames <---
%  set(ud.SpecFilename,'String','C:\My Documents\Boo Boo\Projects and Work\SLIM and GSLIM\SLIM Leg T2\March 10, 2000\Jeff_Te035ms_Chunk002.hdf');
%  set(ud.MaskFilename,'String','C:\My Documents\Boo Boo\Projects and Work\SLIM and GSLIM\SLIM Leg T2\March 10, 2000\CompartmentImage.hdf');
%  set(ud.VarFilename,'String','C:\My Documents\Boo Boo\Projects and Work\SLIM and GSLIM\SLIM Leg T2\March 10, 2000\Jeff_Te035ms.hdf');
%  set(ud.SpecFilename,'String','C:\My Documents\Jeff\Projects and Work\SLIM and GSLIM\SLIM Leg T2\March 10, 2000\Jeff_Te035ms_Chunk002.hdf');
%  set(ud.MaskFilename,'String','C:\My Documents\Jeff\Projects and Work\SLIM and GSLIM\SLIM Leg T2\March 10, 2000\CompartmentImage.hdf');
%  set(ud.VarFilename,'String','C:\My Documents\Jeff\Projects and Work\SLIM and GSLIM\SLIM Leg T2\March 10, 2000\Jeff_Te035ms.hdf');

  % Fill in filename for 2D data
  %set(ud.SpecFilename,'String','/users/bmrl/jtsao/data/GSLIM/Kmiecik/GSLIM/MadeSense/CsiData_dqc3d331_DCinCenter_Slice5_chunk8x8.hdf');
  %set(ud.MaskFilename,'String','/users/bmrl/jtsao/data/GSLIM/Kmiecik/GSLIM/MadeSense/MaskSlice5_fromJoe.hdf');
  
  % Fill in filename for 3D data
  %set(ud.SpecFilename,'String','/users/bmrl/jtsao/data/GSLIM/Kmiecik/GSLIM/MadeSense/CsiData_dqc3d331_DCinCenter_chunk8x8x8.hdf');
  %set(ud.MaskFilename,'String','/users/bmrl/jtsao/data/GSLIM/Kmiecik/GSLIM/MadeSense/Mask.hdf');
  %------------------------------------------------------------------
  
  %------------------------------------------------------------------
  % SET UP DATA
  switch ndims(data),
  case 2,
    if size(data,1)==1,
      data=permute(data,[1,3,4,2]);
    elseif size(data,2)==1,
      data=permute(data,[2,3,4,1]);
    end
  case 3, % Assume last dimension is spectroscopic
    data=permute(data,[1,2,4,3]);
  end
  ud.CsiData = data;
  clear data
  ud.PosX = bitshift(size(ud.CsiData,1),-1)+1;
  ud.PosY = bitshift(size(ud.CsiData,2),-1)+1;
  ud.PosZ = bitshift(size(ud.CsiData,3),-1)+1;
  ud.PosF = bitshift(size(ud.CsiData,4),-1)+1;
  ud.PosXY_h = [];
  ud.PosF_h = [];
  ud.FreqOffset = 4.7;
  ud.SpecWidth = 10;
  ud.ImageAxes = [];
  ud.Image_h = [];
  ud.GraphAxes = [];
  ud.Markers = {};
  set(ud.InputSpecWidth_h,'String',ud.SpecWidth);
  set(ud.InputFreqOffset_h,'String',ud.FreqOffset);
  if isempty(dataname), ud.varname = 'data'; else ud.varname = dataname; end
  %------------------------------------------------------------------
              
  %------------------------------------------------------------------
  % SET UP PAGE
  set(MainFig,'UserData',ud); % Store updated data
  clear ud
  %------------------------------------------------------------------

  %------------------------------------------------------------------
  % Make window visible
  %  SetUpPage(MainFig);		% set up page according to ud.PageNumNow
  SetUpImage(MainFig);
  SetUpGraph(MainFig)
  set(MainFig, 'visible','on','HandleVisibility','callback');
  %------------------------------------------------------------------

  clear Std Btn Edt FigColor MainFig data ud
return
%####################################################################
%####################################################################
%####################################################################
%------- IMAGE CONTROL FUNCTIONS ------------------------------------
function SetUpImage(MainFig)
  if nargin<1,
    MainFig = gcbf;
    if isempty(MainFig), MainFig = gcf; end
  end
  ud = get(MainFig, 'UserData');
  
  %
  % Set up image
  %
  axisH  = ud.ImageAxes;
  imageH = ud.Image_h;
  if isempty(axisH),   % If axis doesn't exist, create it.
    axisH = axes('Interruptible', 'off',...
                   'BusyAction', 'queue',...
                   'Units', 'Pixels',...
                   'Parent', MainFig,...
                   'CLim', [0 1],...
                   'XTick', [],...
                   'YTick', [],...
                   'Position', [2 345-256 256 256],...
                   'Tag', '');
    ud.ImageAxes = axisH;
    set(MainFig,'UserData',ud); % Store updated data
  end
  set(axisH,'ydir', 'reverse', 'XLim',[0 129], 'YLim',[0 129]);
  if isempty(imageH),
    delete(get(axisH,'Children'));    % Delete all children of the axes
    imageH = image('Interruptible', 'off', ...   % Image doesn't exist, creat it.
          'BusyAction', 'queue', ...
          'CData', [], ...
          'Xdata', [1 128], ...
          'Ydata', [1 128], ...
          'CDataMapping', 'Scaled', ...
          'Erasemode', 'none', ...
          'Tag', '', ...
          'Parent', axisH);
    ud.Image_h = imageH;
  else
    childrenH = get(axisH,'Children');
    delete(childrenH(find(childrenH~=imageH)));  % Delete all children except for image
    clear childrenH;
  end
  axes(axisH);  % set as current axis  
  title(''); xlabel('');
  set(MainFig,'UserData',ud); % Store updated data
%clear axisH imageH         % If a new image has been created, imageH is no long up-to-date.

  if isempty(ud.CsiData) | ud.PosF<1 | ud.PosF>size(ud.CsiData,4),
    ud.DisplayDataSet = 0;
    set(MainFig,'UserData',ud); % Store updated data
    % Update image control
    set(ud.SliceNum_h, 'String', '');
    set(ud.ImageControls_h, 'Enable', 'off');
    set(ud.ImageSliceCtrls_h, 'Enable', 'off');
    
    % Update image display
    set(ud.Image_h,'CData',0);
    caxis([-1, 1]);  % set display at intermediate color. Otherwise, looks like there is an image, but it's too dim to see.
  else
    ud.DisplayDataSet = double(ud.CsiData(:,:,:,ud.PosF));	% Summing along spectroscopic dimension
    % Pad X and Y dimension to at least 33 for better visualization
    % if (size(ud.DisplayDataSet,1)<33 | size(ud.DisplayDataSet,2)<33),
      % Zeropad in transform space,
      % This performs the same function as Jeff's ShowSpecData_DCinCenter
      SizeX = size(ud.DisplayDataSet,1);
      SizeY = size(ud.DisplayDataSet,2);
      SizeZ = size(ud.DisplayDataSet,3);
      
      if SizeX > 1,
        set(axisH,'XLim',[0.5 SizeX+0.5]);
        set(imageH,'Xdata',[0.5 SizeX+0.5]);
      else
        set(axisH,'XLim',[0.5 1.5]);
        set(imageH,'Xdata',[0.5 1.5]);
      end
      if SizeY > 1,
        set(axisH,'YLim',[0.5 SizeY+0.5]);
        set(imageH,'Ydata',[0.5 SizeY+0.5]);
      else
        set(axisH,'YLim',[0.5 1.5]);
        set(imageH,'Ydata',[0.5 1.5]);
      end
 
      % Zeropadding ...
      ZeroPadMinSize = 65;
      if SizeX<ZeroPadMinSize & SizeX>1,
        % Equivalent to shift DC to corner -> IFFT, shift DC back to center
          p = bitshift(SizeX,-1)+1;
          ud.DisplayDataSet([1:SizeX],:,:) = ud.DisplayDataSet([p:SizeX 1:p-1],:,:);	% PI Shift in X
          ud.DisplayDataSet = ifft(ud.DisplayDataSet,[],1);	% ifft2 is broken, so fill use 1D FT all the time
          ud.DisplayDataSet([p:SizeX 1:p-1],:,:) = ud.DisplayDataSet([1:SizeX],:,:);	% PI Shift in K
          % Zero pad
          tmpDisplayDataSet = zeros(ZeroPadMinSize,SizeY,SizeZ);
          tmpDisplayDataSet([1:SizeX]-bitshift(SizeX,-1)+bitshift(end,-1),:,:) = ud.DisplayDataSet;
          ud.DisplayDataSet = tmpDisplayDataSet; clear tmpDisplayDataSet
          SizeX = ZeroPadMinSize;
          
          % Equivalent to shift DC to corner -> FFT, shift DC back to center
          p = bitshift(SizeX,-1)+1;
          ud.DisplayDataSet([1:SizeX],:,:) = ud.DisplayDataSet([p:SizeX 1:p-1],:,:);	% PI Shift in X
          ud.DisplayDataSet = fft(ud.DisplayDataSet,[],1);	% ifft2 is broken, so fill use 1D FT all the time
          ud.DisplayDataSet([p:SizeX 1:p-1],:,:) = ud.DisplayDataSet([1:SizeX],:,:);	% PI Shift in Kx
          clear p
      end
      if SizeY<ZeroPadMinSize & SizeY>1,
          % PI shift -> FFT -> PI shift
          % Equivalent to shift DC to corner -> FFT, shift DC back to center
          p = bitshift(SizeY,-1)+1;
          ud.DisplayDataSet(:,[1:SizeY],:) = ud.DisplayDataSet(:,[p:SizeY 1:p-1],:);	% PI Shift in Ky
          ud.DisplayDataSet = ifft(ud.DisplayDataSet,[],2);
          ud.DisplayDataSet(:,[p:SizeY 1:p-1],:) = ud.DisplayDataSet(:,[1:SizeY],:);	% PI Shift in Y
          
  		    % Zero pad
          tmpDisplayDataSet = zeros(SizeX,ZeroPadMinSize,SizeZ);
          tmpDisplayDataSet(:,[1:SizeY]-bitshift(SizeY,-1)+bitshift(end,-1),:) = ud.DisplayDataSet;
          ud.DisplayDataSet = tmpDisplayDataSet; clear tmpDisplayDataSet
          SizeY = ZeroPadMinSize;

          % Equivalent to shift DC to corner -> IFFT, shift DC back to center
          p = bitshift(SizeY,-1)+1;
          ud.DisplayDataSet(:,[1:SizeY],:) = ud.DisplayDataSet(:,[p:SizeY 1:p-1],:);	% PI Shift in Y
          ud.DisplayDataSet = fft(ud.DisplayDataSet,[],2);
          ud.DisplayDataSet(:,[p:SizeY 1:p-1],:) = ud.DisplayDataSet(:,[1:SizeY],:);	% PI Shift in Ky
          clear p
        end       
      clear SizeX SizeY SizeZ
    % end	 % Finished zeropadding
    clear ZeroPadMinSize    
    %% Take absolute value
    %ud.DisplayDataSet = abs(ud.DisplayDataSet);   
    
    %ud.PosX = bitshift(size(dataset,1),-1)+1;
    %ud.PosY = bitshift(size(dataset,2),-1)+1;
    set(MainFig,'UserData',ud); % Store updated data
  
    % Update image control
    set(ud.ImageControls_h, 'Enable', 'on');
    if size(ud.CsiData,3)>1,
      % If there is >1 slice, enable image slice control
      set(ud.ImageSliceCtrls_h, 'Enable', 'on');
    else
      set(ud.ImageSliceCtrls_h, 'Enable', 'off');
    end

    % Update image display
    ShowImageSlice(ud.PosZ, MainFig); %bitshift(size(ud.DisplayDataSet,3),-1)+1, MainFig);
    SetUpColorScale(0, MainFig);   % Set up color scale
    SetUpSpatPos(MainFig);
  end % end of choosing between no data or have data
  clear ud
  clear axisH imageH
return
%####################################################################
%####################################################################
%####################################################################
function SetUpColorScale(refreshimage, MainFig)
  if nargin<1, refreshimage = 0; end
  if nargin<2,
    MainFig = gcbf;
    if isempty(MainFig), MainFig = gcf; end
  end
  ud = get(MainFig, 'UserData');

  % Find minimum and maximum of image intensity
  whichVal = get(ud.ImageComplex_h,'value');
  if whichVal == 4,     % Phase
    minval=-pi; maxval=pi;
  else
    switch whichVal,
      case 2,    tempImg = real(ud.DisplayDataSet(:));
      case 3,    tempImg = imag(ud.DisplayDataSet(:));
      otherwise, tempImg =  abs(ud.DisplayDataSet(:));
    end
    maxval = max(tempImg);
    minval = min(tempImg);
    clear tempImg
    if maxval<=minval, maxval=minval+1; minval=minval-1; end % make sure that intensity range is valid
  end
  clear whichVal
  
  % Set image axes as current axes
  if ~isempty(ud.ImageAxes), axes(ud.ImageAxes); end % set as current axis
  
  % Update colorbar
  caxisNow = caxis;
  if caxisNow(1)~=minval | caxisNow(2)~=maxval, caxis([minval, maxval]); end
  clear caxisNow minval maxval
  
  % Refresh image
  if refreshimage, GoToSlice('refresh', MainFig); end
  clear ud
return
%####################################################################
%####################################################################
%####################################################################
function ShowImageSlice(slidenum, MainFig)
  if nargin<1,
    return
  end
  if nargin<2,
    MainFig = gcbf;
    if isempty(MainFig), MainFig = gcf; end
  end
  
  ud = get(MainFig, 'UserData');
  
  % Update image
  switch get(ud.ImageComplex_h,'value'),
    case 2
      tempImg = real(ud.DisplayDataSet(:,:,slidenum));
    case 3
      tempImg = imag(ud.DisplayDataSet(:,:,slidenum));
    case 4
      tempImg = angle(ud.DisplayDataSet(:,:,slidenum));
    otherwise
      tempImg = abs(ud.DisplayDataSet(:,:,slidenum));
  end
  set(ud.Image_h,'CData',transpose(tempImg));
  % zoom on
  clear tempImg
  
  % Update slice number display
  set(ud.SliceNum_h, 'String', ...
     {sprintf('Slice %d (1 to %d)',slidenum,size(ud.DisplayDataSet,3))});
  ud.PosZ = slidenum;
  
  set(MainFig,'UserData',ud); % Store updated data
  clear ud
return
%####################################################################
%####################################################################
%####################################################################
function GoToSlice(goWhere, MainFig)
  if nargin<1,
    return
  end
  if nargin<2,
    MainFig = gcbf;
    if isempty(MainFig), MainFig = gcf; end
  end

  ud = get(MainFig, 'UserData');
  slidenum = ud.PosZ;
  LastSlide = size(ud.DisplayDataSet,3);
  switch goWhere
    case 'first'
      if slidenum == 1, return; end
      ud.PosZ = 1;
    case 'prev'
      if slidenum == 1, return; end
      ud.PosZ = ud.PosZ-1;
    case 'next'
      if slidenum == LastSlide, return; end
      ud.PosZ = ud.PosZ+1;
    case 'last'
      if slidenum == LastSlide, return; end
      ud.PosZ = LastSlide;
    case 'refresh' % keep the same PosZ
   end
  ShowImageSlice(ud.PosZ, MainFig);
  SetUpSpatPos(MainFig);
  SetUpGraph(MainFig);
  clear ud slidenum LastSlide
return

%####################################################################
%####################################################################
%####################################################################
%------- COLOR PALETTE CONTROL FUNCTIONS ----------------------------
function changeColors(colorname)
  if nargin<1,
    return
  end
  try, colormap(lower(colorname)); end
return
%####################################################################
%####################################################################
%####################################################################
function invertColors()
  colormap(flipud(colormap)); 
return
%####################################################################
%####################################################################
%####################################################################
function changeBrightness(amount)
  if nargin<1,
    return
  end
  brighten(amount);
return
%####################################################################
%####################################################################
%####################################################################
function SetUpGraph(MainFig)
  if nargin<1, MainFig = gcbf; if isempty(MainFig), MainFig = gcf; end; end
  ud = get(MainFig, 'UserData');
  %
  % Set up image
  %
  axisH  = ud.GraphAxes;
  if isempty(axisH),   % If axis doesn't exist, create it.
    axisval='';
    axisH = axes('Interruptible', 'off',...
                   'BusyAction', 'queue',...
                   'Units', 'Pixels',...
                   'Parent', MainFig,...
                   'ydir', 'reverse',...
                    'CLim', [0 1],...
                   'XTick', [],...
                   'YTick', [],...
                   'Position', [2+300 345-256 256 256],...
                   'Tag', '');
    axes(axisH);
    ud.GraphAxes = axisH;
    set(MainFig,'UserData',ud); % Store updated data
  else
    axes(axisH);  % set as current axis
    axisval = axis;
    set(axisH,'ydir', 'reverse');
    delete(get(axisH,'Children'));    % Delete all children
  end
  
  if isempty(ud.CsiData),
    plot(0,0);
  else
    len = size(ud.CsiData,4);
    x=([1:len]-(bitshift(len,-1)+1))./len*ud.SpecWidth + ud.FreqOffset;
    hold off;
    if get(ud.PlotImag_h,'Value'),
       plot(x,reshape(imag(double(ud.CsiData(ud.PosX,ud.PosY,ud.PosZ,:))), len, 1), 'Color', [1,0,0]);
       hold on;
    end
    if get(ud.PlotReal_h,'Value'),
       plot(x,reshape(real(double(ud.CsiData(ud.PosX,ud.PosY,ud.PosZ,:))), len, 1), 'Color', [0,0,1]);
       hold on;
    end
    if get(ud.PlotAbs_h,'Value'),
       plot(x,reshape(abs(double(ud.CsiData(ud.PosX,ud.PosY,ud.PosZ,:))), len, 1), 'Color', [0,0,0]);
    end
    hold off;
    set(axisH, 'xdir', 'reverse');
    axis tight;
    if ~isempty(axisval),
       newaxisval=axis;
       axis([newaxisval([1,2]),axisval([3,4])]);
       clear newaxisval;
    end
    SetUpSpecPos(MainFig);
    clear len h
  end
  clear ud axisval axisH
  xlabel('');
return
%####################################################################
%####################################################################
%####################################################################
function SetUpSpecPos(MainFig)
  if nargin<1, MainFig = gcbf; if isempty(MainFig), MainFig = gcf; end; end
  ud = get(MainFig, 'UserData');
  if ~isempty(ud.PosF_h), try, delete(ud.PosF_h); end; end
  
  axisH = ud.GraphAxes;
  if ~isempty(axisH),
    axes(axisH);
    x=(ud.PosF-(bitshift(size(ud.CsiData,4),-1)+1))./size(ud.CsiData,4)*ud.SpecWidth + ud.FreqOffset;
    ud.PosF_h = line([x,x],get(axisH,'YLim'),'Color', [0,0.6,0],'LineStyle',':','Tag','');
    clear x
    title(sprintf('Pos %d (1 to %d): %.2fppm  at (%d,%d,%d)',...
       ud.PosF,size(ud.CsiData,4),(ud.PosF-(bitshift(size(ud.CsiData,4),-1)+1))/size(ud.CsiData,4)*ud.SpecWidth + ud.FreqOffset,ud.PosX,ud.PosY,ud.PosZ));
    set(MainFig,'UserData',ud); % Store updated data
  end
return
%####################################################################
%####################################################################
%####################################################################
function SetUpSpatPos(MainFig)
  if nargin<1, MainFig = gcbf; if isempty(MainFig), MainFig = gcf; end; end
  ud = get(MainFig, 'UserData');
  
  if ~isempty(ud.PosXY_h), try, delete(ud.PosXY_h); end; end
  
  if ~isempty(ud.ImageAxes),
    axes(ud.ImageAxes);
    ud.PosXY_h = line(ud.PosX,ud.PosY,'Parent',ud.ImageAxes,'Color',[0,0.8,0],...
       'Marker','+','MarkerSize',15, 'Tag','');
    set(MainFig,'UserData',ud); % Store updated data
  end
return
%####################################################################
%####################################################################
%####################################################################
function TrackMouseMotion(MainFig)
  if nargin<1,
    MainFig = gcbf;
    if isempty(MainFig), MainFig = gcf; end
  end
  ud = get(MainFig, 'UserData');
  pt = get(MainFig, 'CurrentPoint');
  x = pt(1,1);
  y = pt(1,2);
  
  pos = get(ud.ImageAxes, 'Position');
  if (x>pos(1) & x<pos(1)+pos(3) & y>pos(2) & y<pos(2)+pos(4))
    set(MainFig, 'Pointer', 'cross');
  else
    pos = get(ud.GraphAxes, 'Position');
    if (x>pos(1) & x<pos(1)+pos(3) & y>pos(2) & y<pos(2)+pos(4))
      set(MainFig, 'Pointer', 'cross');
    else
      set(MainFig, 'Pointer', 'arrow');
    end
  end
  clear x y pt
return
%####################################################################
%####################################################################
%####################################################################
function TrackMouseButton(MainFig)
  if nargin<1,
    MainFig = gcbf;
    if isempty(MainFig), MainFig = gcf; end
  end
  ud = get(MainFig, 'UserData');
  pt = get(MainFig, 'CurrentPoint');
  x = pt(1,1);
  y = pt(1,2);
  
  axisH = ud.ImageAxes;
  pos = get(axisH, 'Position');
  if (x>pos(1) & x<pos(1)+pos(3) & y>pos(2) & y<pos(2)+pos(4)) % In left image
    pt = get(axisH, 'CurrentPoint');
    x = round(pt(1,1));
    y = round(pt(1,2));
    if     x<1,                  ud.PosX = 1;
    elseif x>size(ud.CsiData,1), ud.PosX = size(ud.CsiData,1);
    else                         ud.PosX = x; end
    if     y<1,                  ud.PosY = 1;
    elseif y>size(ud.CsiData,2), ud.PosY = size(ud.CsiData,2);
    else                         ud.PosY = y; end
    set(MainFig,'UserData',ud); % Store updated data
    SetUpGraph(MainFig);
    SetUpSpatPos(MainFig);
    clear axisH pos x y pt ud
    return
  end
  
  axisH = ud.GraphAxes;
  pos = get(axisH, 'Position');
  if (x>pos(1) & x<pos(1)+pos(3) & y>pos(2) & y<pos(2)+pos(4))
    pt = get(axisH, 'CurrentPoint');
    x = (pt(1,1)-ud.FreqOffset)/ud.SpecWidth*size(ud.CsiData,4)+(bitshift(size(ud.CsiData,4),-1)+1);
    x = round(x);
    if     x<1, x=1;
    elseif x>size(ud.CsiData,4), x=size(ud.CsiData,4);
    end
    if ud.PosF ~= x,
      ud.PosF = x;
      set(MainFig,'UserData',ud); % Store updated data
      SetUpImage(MainFig);
      SetUpSpecPos(MainFig);
    end
    clear axisH pos x y pt ud
    return
  end
return
%####################################################################
%####################################################################
%###################################################################
function ChangeFreqPos(ChangeAmount, MainFig)
  if nargin<1, return; end % No change given
  if nargin<2,
    MainFig = gcbf;
    if isempty(MainFig), MainFig = gcf; end
  end
  ud = get(MainFig, 'UserData');
  if size(ud.CsiData,4)<=1, return; end % Only 1 frequency, no change possible
  
  if ud.SpecWidth<0, ChangeAmount = -ChangeAmount; end  % Change direction of change according to sign of spectral bandwidth
  
  % Change the frequency position
  if ChangeAmount<0,
    if ud.PosF==1, return; end % Already reached the start
    if ud.PosF+ChangeAmount<1, ud.PosF = 1;
    else ud.PosF = ud.PosF+ChangeAmount; end
  else
    if  ud.PosF==size(ud.CsiData,4), return; end % Already reached the end
    if ud.PosF+ChangeAmount>size(ud.CsiData,4),
       ud.PosF = size(ud.CsiData,4);
    else ud.PosF = ud.PosF+ChangeAmount; end
  end
  
  % Update image and frequency position
  set(MainFig,'UserData',ud); % Store updated data
  SetUpImage(MainFig);
  SetUpSpecPos(MainFig);
return  
%####################################################################
%####################################################################
%###################################################################
function changeSpecWidth(newval, MainFig)
  if nargin<1, return;  end
  if nargin<2,
    MainFig = gcbf;
    if isempty(MainFig), MainFig = gcf; end
  end
  ud = get(MainFig, 'UserData');
  newval = eval(newval, num2str(ud.SpecWidth));
  if ~isempty(ud.InputSpecWidth_h), set(ud.InputSpecWidth_h, 'String', num2str(newval), 'UserData', newval); end
  
  ud.SpecWidth = newval;
  set(MainFig,'UserData',ud); % Store updated data
  clear ud
  SetUpGraph(MainFig);
return
%####################################################################
%####################################################################
%####################################################################
function changeFreqOffset(newval, MainFig)
  if nargin<1, return;  end
  if nargin<2,
    MainFig = gcbf;
    if isempty(MainFig), MainFig = gcf; end
  end
  ud = get(MainFig, 'UserData');
  newval = eval(newval, num2str(ud.SpecWidth));
  if ~isempty(ud.InputFreqOffset_h), set(ud.InputFreqOffset_h, 'String', num2str(newval), 'UserData', newval); end
  
  ud.FreqOffset = newval;
  set(MainFig,'UserData',ud); % Store updated data
  clear ud
  SetUpGraph(MainFig);
return
%####################################################################
%####################################################################
%####################################################################
function ProcessData(action, MainFig)
  if nargin<1, return;  end
  if nargin<2,
    MainFig = gcbf;
    if isempty(MainFig), MainFig = gcf; end
  end
  ud = get(MainFig, 'UserData');
  ApplyToWhichDim = [get(ud.ApplyToDimXh,'Value'),...
                     get(ud.ApplyToDimYh,'Value'),...
                     get(ud.ApplyToDimZh,'Value'),...
                     get(ud.ApplyToDimFh,'Value')];
  if isempty(find(ApplyToWhichDim~=0)), return; end  % No dimension chosen                 
  
  switch action,
    case 'FFT' %---------------------------------------
      idx = find(ApplyToWhichDim~=0);
      ud.CsiData = transformImageToKspace(ud.CsiData,idx);
      clear idx
    case 'IFFT' %---------------------------------------
      idx = find(ApplyToWhichDim~=0);
      ud.CsiData = transformKspaceToImage(ud.CsiData,idx);
      clear idx
    case 'Shift to center' %---------------------------------------
      ShiftAmount = ( bitshift([size(ud.CsiData,1),size(ud.CsiData,2),size(ud.CsiData,3),size(ud.CsiData,4)],-1)+1 ...
                    - [ud.PosX,ud.PosY,ud.PosZ,ud.PosF]) .*(ApplyToWhichDim~=0);
      if sum(abs(ShiftAmount))==0, clear ShiftAmount; return; end
      ud.CsiData = shift(ud.CsiData,ShiftAmount);
      if ApplyToWhichDim(1), ud.PosX=ud.PosX+ShiftAmount(1); end
      if ApplyToWhichDim(2), ud.PosY=ud.PosY+ShiftAmount(2); end
      if ApplyToWhichDim(3), ud.PosZ=ud.PosZ+ShiftAmount(3); end
      if ApplyToWhichDim(4), ud.PosF=ud.PosF+ShiftAmount(4); end
    case 'Zeropad to center' %---------------------------------------
      CurrentPos = [ud.PosX,ud.PosY,ud.PosZ,ud.PosF];
      OrigSize = {size(ud.CsiData,1),size(ud.CsiData,2),size(ud.CsiData,3),size(ud.CsiData,4)};
      NewSize = OrigSize;
      ShiftAmount = zeros(1,4);
      for n=[1:4],
        if ApplyToWhichDim(n),
          if bitshift(size(ud.CsiData,n),-1)+1<CurrentPos(n), % Pad at the end
            NewSize{n} = 2*CurrentPos(n)-2;
          elseif bitshift(size(ud.CsiData,n),-1)+1>CurrentPos(n), % Pad in front
            NewSize{n} = 2*(size(ud.CsiData,n) - CurrentPos(n)) + 1;
            ShiftAmount(n) = NewSize{n}-size(ud.CsiData,n);
          end
        end
      end
      if isequal(NewSize,OrigSize), return; end
      ud.CsiData(NewSize{:}) = 0; % Zero-pad
      if ~isequal(ShiftAmount, zeros(1,4)),
        ud.CsiData = shift(ud.CsiData, ShiftAmount); % Apply circular shift
      end
      if ApplyToWhichDim(1), ud.PosX=ud.PosX+ShiftAmount(1); end
      if ApplyToWhichDim(2), ud.PosY=ud.PosY+ShiftAmount(2); end
      if ApplyToWhichDim(3), ud.PosZ=ud.PosZ+ShiftAmount(3); end
      if ApplyToWhichDim(4), ud.PosF=ud.PosF+ShiftAmount(4); end
      clear n NewSize OrigSize ShiftAmount
    case 'Multiply with Hamming' %---------------------------------------
      for n=[1:4], idx{n} = ones(1,size(ud.CsiData,n)); end
      for n=[1:4],
         if ApplyToWhichDim(n),
          filter = shiftdim( sqrt(hamming(size(ud.CsiData,n))), 1-n);
          idx{n} = [1:size(ud.CsiData,n)];
          ud.CsiData = ud.CsiData.*filter(idx{:});
          idx{n} = ones(1,size(ud.CsiData,n));
          clear filter
        end
      end
      clear n idx
    case 'Divide by Hamming' %---------------------------------------
      for n=[1:4], idx{n} = ones(1,size(ud.CsiData,n)); end
      for n=[1:4],
        if ApplyToWhichDim(n),
          filter = shiftdim( 1./sqrt(hamming(size(ud.CsiData,n))), 1-n);
          idx{n} = [1:size(ud.CsiData,n)];
          ud.CsiData = ud.CsiData.*filter(idx{:});
          idx{n} = ones(1,size(ud.CsiData,n));
          clear filter
        end
      end
      clear n idx
    case 'Crop Positive' %---------------------------------------
      CurrentPos = [ud.PosX,ud.PosY,ud.PosZ,ud.PosF];
      OrigSize = [size(ud.CsiData,1),size(ud.CsiData,2),size(ud.CsiData,3),size(ud.CsiData,4)];
      idx = cell(1, 4);
      AnythingChanged = 0;
      for n=[1:4],
        if ApplyToWhichDim(n), 
          if n==4 & ud.SpecWidth<0,
            if CurrentPos(n)~=1, AnythingChanged = 1; end
            idx{n} = [CurrentPos(n):OrigSize(n)];
          else
            if CurrentPos(n)~=OrigSize(n), AnythingChanged = 1; end
            idx{n} = [1:CurrentPos(n)];
          end
        else
          idx{n} = [1:OrigSize(n)];
        end
      end; clear n
      if AnythingChanged == 0, return; end % Nothing to change
      ud.CsiData = ud.CsiData(idx{:});
      if ApplyToWhichDim(4) & ud.SpecWidth<0, ud.PosF = 1; end
      FreqPos = ([min(idx{4}), max(idx{4})]-(bitshift(OrigSize(4),-1)+1))./OrigSize(4)*ud.SpecWidth + ud.FreqOffset;
      ud.SpecWidth = (FreqPos(2)-FreqPos(1))/(size(ud.CsiData,4)-1)*size(ud.CsiData,4);
      ud.FreqOffset = FreqPos(1) + bitshift(size(ud.CsiData,4),-1)./size(ud.CsiData,4)*ud.SpecWidth;
      set(ud.InputSpecWidth_h,'String',ud.SpecWidth);
      set(ud.InputFreqOffset_h,'String',ud.FreqOffset);
      clear CurrentPos OrigSize n AnythingChanged idx FreqPos
    case 'Crop Negative' %---------------------------------------
      CurrentPos = [ud.PosX,ud.PosY,ud.PosZ,ud.PosF];
      OrigSize = [size(ud.CsiData,1),size(ud.CsiData,2),size(ud.CsiData,3),size(ud.CsiData,4)];
      idx = cell(1, 4);
      AnythingChanged = 0;
      for n=[1:4],
        if ApplyToWhichDim(n), 
          if n==4 & ud.SpecWidth<0,
            if CurrentPos(n)~=OrigSize(n), AnythingChanged = 1; end
            idx{n} = [1:CurrentPos(n)];
          else
            if CurrentPos(n)~=1, AnythingChanged = 1; end
            idx{n} = [CurrentPos(n):OrigSize(n)];
          end
        else
          idx{n} = [1:OrigSize(n)];
        end
      end; clear n
      if AnythingChanged == 0, return; end % Nothing to change
      ud.CsiData = ud.CsiData(idx{:});
      if ApplyToWhichDim(1), ud.PosX = 1; end
      if ApplyToWhichDim(2), ud.PosY = 1; end
      if ApplyToWhichDim(3), ud.PosZ = 1; end
      if ApplyToWhichDim(4) & ud.SpecWidth>0, ud.PosF = 1; end
      FreqPos = ([min(idx{4}), max(idx{4})]-(bitshift(OrigSize(4),-1)+1))./OrigSize(4)*ud.SpecWidth + ud.FreqOffset;
      ud.SpecWidth = (FreqPos(2)-FreqPos(1))/(size(ud.CsiData,4)-1)*size(ud.CsiData,4);
      ud.FreqOffset = FreqPos(1) + bitshift(size(ud.CsiData,4),-1)./size(ud.CsiData,4)*ud.SpecWidth;
      set(ud.InputSpecWidth_h,'String',ud.SpecWidth);
      set(ud.InputFreqOffset_h,'String',ud.FreqOffset);
      clear CurrentPos OrigSize n AnythingChanged FreqPos idx
    otherwise
      fprintf('Unknown action in ProcessData: %s\n', action);
      return
  end
    
  set(MainFig,'UserData',ud); % Store updated data
  SetUpImage(MainFig);
  SetUpGraph(MainFig);
    
    %  newval = eval(newval, num2str(ud.SpecWidth));
%  h=ud.InputSpecWidth_h;
%  if ~isempty(h), set(h, 'String', num2str(newval), 'UserData', newval); end
%  
 % ud.SpecWidth = newval;
%  set(MainFig,'UserData',ud); % Store updated data
%  clear ud
%  SetUpGraph(MainFig);
return
%####################################################################
%####################################################################
%####################################################################
function AutoScaleAxis(whichAxis, MainFig)
  if nargin<1, return; end
  if nargin<2, MainFig = gcbf; if isempty(MainFig), MainFig = gcf; end; end
  ud = get(MainFig, 'UserData');
  
  axisH  = ud.GraphAxes;
  if ~isempty(find(whichAxis=='x')), set(axisH,'XLimMode','auto'); end
  if ~isempty(find(whichAxis=='y')), set(axisH,'YLimMode','auto'); end
  SetUpSpecPos(MainFig);
return
%####################################################################
%####################################################################
%####################################################################
function AskForVariableName(MainFig)
  if nargin<1, MainFig = gcbf; if isempty(MainFig), MainFig = gcf; end; end
  ud = get(MainFig, 'UserData');
  bmf = get(MainFig,'WindowButtonMotionFcn');
  bdf = get(MainFig,'WindowButtonDownFcn');
  set(MainFig,'WindowButtonMotionFcn','','WindowButtonDownFcn','');
  
  h = dialog('Visible','off','WindowStyle','modal');

  set(h,'Units','pixels','Position',[100 300 320 70],'Visible','on','Tag','polydlg');
  uicontrol(h,'Style','text','Units','pixels','Position',[10 40 200 20],...
          'Horiz','left','String','Name of variable to save data to:');

  ud.varnameH = uicontrol(h,'Style','edit','Units','pixels',...
     'String',ud.varname,'HorizontalAlignment','left','Position',[210 40 100 20],...
     'BackgroundColor', [1,1,1]);
  uicontrol('Style','Pushbutton','Units','pixels','Position',[200 10 50 20],'String','OK',...
     'Callback', 'MainFig = get(gco,''Parent''); ud=get(MainFig,''UserData''); varname=get(ud.varnameH,''string''); if ~isempty(varname), ud.varname=varname; set(MainFig,''UserData'', ud); assignin(''base'',varname,ud.CsiData); end; clear MainFig ud varname; close(gcbf)');
  uicontrol('Style','Pushbutton','Units','pixels','Position',[260 10 50 20],...
     'Callback','close(gcbf)','String','Cancel');
  set(h,'HandleVisibility','off')
  set(MainFig,'UserData',ud,'WindowButtonMotionFcn',bmf,'WindowButtonDownFcn',bdf);
return
%####################################################################
%####################################################################
%####################################################################
function w = hamming(n,sflag)
%HAMMING Hamming window.
%   W = HAMMING(N) returns the N-point symmetric Hamming window 
%       in a column vector. 
%   W = HAMMING(N,SFLAG) generates the N-point Hamming window 
%       using SFLAG window sampling. SFLAG may be either 'symmetric' 
%       or 'periodic'. By default, 'symmetric' window sampling is used. 
error(nargchk(1,2,nargin));
if n==1, w=1; return; end

% Set sflag to default if it's not already set:
if nargin == 1, sflag = 'symmetric'; end

switch lower(sflag),
case 'periodic'
   w = .54 - .46*cos(2*pi*(0:n-1)'/(n));
case 'symmetric'
   w = .54 - .46*cos(2*pi*(0:n-1)'/(n-1));
otherwise
	error('Sampling must be either ''symmetric'' or ''periodic''');
end
%####################################################################
%####################################################################
%####################################################################
% VERSION:
%   1.0 - March 27, 2000 Jeffrey Tsao
%
% [newdata] = transformImageToKspace(data,dim,DCasFirstElement,NormalizeAsInV)
%
% - To transform an image to k-space by Forward Fourier Transform.
% - Input Parameters
%     data: image
%     dim: dimension(s) to transform, if not given, transform all dimensions
%     DCasFirstElement: by default, the DC term of image and k-space are 
%                       considered to be in the center (bitshift(datasize,-1)+1).
%                       If DCasFirstElement==1, DC is considered to be at the first
%                       element.
%     NormalizeAsInV: 1 to normalize Fourier transform according to the
%                     convention of V (default). 0 to use native convention, which
%                     may be variable from computer to computer.
%
function [newdata] = transformImageToKspace(data,dim,DCasFirstElement,NormalizeAsInV)
  if nargin<1, help mfilename; end
  if nargin<2, dim=[]; end
  if nargin<3, DCasFirstElement=0; end
  if nargin<4, NormalizeAsInV=1; end
  
  numDims = ndims(data);   % Number of dimensions

  %----------------------------------------------------------------------
  % Move DC from center to 1st element, if needed.
  %   (fft routine assumes DC to be 1st element)
  if DCasFirstElement,
    newdata = data;
  else                          % Move DC to 1st element for fft
    idx = cell(1, numDims);
    for k = 1:numDims
      m = size(data, k);
      if m>1 & (isempty(dim) | ~isempty(find(k==dim))), % either transform all dimensions or just transform this dimension
        p = bitshift(m,-1)+1;   % central pixel (i.e. position of DC)
        idx{k} = [p:m 1:p-1];
        clear p
      else
        idx{k} = [1:m];
      end
    end
    clear k m
    newdata = data(idx{:});      % Perform fft-shift
  end
  %----------------------------------------------------------------------
  
  %----------------------------------------------------------------------
  % Perform fft
  fftnumelements = 1;
  for k = 1:numDims
    m = size(newdata, k);
    if m>1 & (isempty(dim) | ~isempty(find(k==dim))), % either transform all dimensions or just transform this dimension
      newdata = fft(newdata,[],k);
      fftnumelements = fftnumelements*m;  % Count total number of elements in fft
    end
  end
  clear k m numDims
  %----------------------------------------------------------------------
  
  %----------------------------------------------------------------------
  % Move DC from 1st element back to center, if necessary
  if DCasFirstElement==0,
    newdata(idx{:}) = newdata;  % Perform fft-shift
    clear idx
  end
  %----------------------------------------------------------------------
  clear numDims
  
  %----------------------------------------------------------------------
  % Normalize according to V convention, if necessary
  if NormalizeAsInV,
    newdata = newdata/fftnumelements;
  end
  %----------------------------------------------------------------------
return
%####################################################################
%####################################################################
%####################################################################
% VERSION:
%   1.0 - March 27, 2000 Jeffrey Tsao
%
% [newdata] = transformKspaceToImage(data,dim,DCasFirstElement,NormalizeAsInV)
%
% - To transform k-space data to an image by Inverse Fourier Transform.
% - Input Parameters
%     data: k-space data
%     dim: dimension(s) to transform, if not given, transform all dimensions
%     DCasFirstElement: by default, the DC term of image and k-space are 
%                       considered to be in the center (bitshift(datasize,-1)+1).
%                       If DCasFirstElement==1, DC is considered to be at the first
%                       element.
%     NormalizeAsInV: 1 to normalize Fourier transform according to the
%                     convention of V (default). 0 to use native convention, which
%                     may be variable from computer to computer.
%
function [newdata] = transformKspaceToImage(data,dim,DCasFirstElement,NormalizeAsInV)
  if nargin<1, help mfilename; end
  if nargin<2, dim=[]; end
  if nargin<3, DCasFirstElement=0; end
  if nargin<4, NormalizeAsInV=1; end
  
  numDims = ndims(data);   % Number of dimensions

  %----------------------------------------------------------------------
  % Move DC from center to 1st element, if needed.
  %   (fft routine assumes DC to be 1st element)
  if DCasFirstElement,
    newdata = data;
  else                          % Move DC to 1st element for fft
    idx = cell(1, numDims);
    for k = 1:numDims
      m = size(data, k);
      if m>1 & (isempty(dim) | ~isempty(find(k==dim))), % either transform all dimensions or just transform this dimension
        p = bitshift(m,-1)+1;   % central pixel (i.e. position of DC)
        idx{k} = [p:m 1:p-1];
        clear p
      else
        idx{k} = [1:m];
      end
    end
    clear k m
    newdata = data(idx{:});      % Perform fft-shift
  end
  %----------------------------------------------------------------------
  
  %----------------------------------------------------------------------
  % Perform fft
  fftnumelements = 1;
  for k = 1:numDims
    m = size(newdata, k);
    if m>1 & (isempty(dim) | ~isempty(find(k==dim))), % either transform all dimensions or just transform this dimension
      newdata = ifft(newdata,[],k);
      fftnumelements = fftnumelements*m;  % Count total number of elements in fft
    end
  end
  clear k m numDims
  %----------------------------------------------------------------------
  
  %----------------------------------------------------------------------
  % Move DC from 1st element back to center, if necessary
  if DCasFirstElement==0,
    newdata(idx{:}) = newdata;  % Perform fft-shift
    clear idx
  end
  %----------------------------------------------------------------------
  clear numDims
  
  %----------------------------------------------------------------------
  % Normalize according to V convention, if necessary
  if NormalizeAsInV,
    newdata = newdata*fftnumelements;
  end
  %----------------------------------------------------------------------
return
%####################################################################
%####################################################################
%####################################################################
%
% [y] = shift(x, offset, useLinearShift)
%
% - Shift data by specified offset
% - Assumed circular shift if "useLinearShift" is not given
% 
function [y] = shift(x, offset, useLinearShift)
  if nargin<3,
    useLinearShift = 0;
  elseif useLinearShift=='y' | useLinearShift=='Y' | useLinearShift==1,
    useLinearShift = 1;
  else
    useLinearShift = 0;
  end
  
  numDims = ndims(x);
  if prod(size(offset))==1,
    offset = ones(1,numDims)*offset;
    offset = offset.*(size(x)~=1);
  end
  
  idx = cell(1, numDims);
  if useLinearShift, recipient = cell(1, numDims);  end
  for k = 1:numDims
    m = size(x, k);
    if k<=prod(size(offset)),
      if useLinearShift,
        if offset(k)>=0,
          idx{k}       = [1:m-offset(k)];
          recipient{k} = [1+offset(k):m];
        else
          idx{k}       = [1-offset(k):m];
          recipient{k} = [1:m+offset(k)];          
        end
      else
        offset(k) = mod(offset(k), size(x,k));
        idx{k} = [m-offset(k)+1:m 1:m-offset(k)];
      end
    else
      idx{k} = [1:m];
      if useLinearShift, recipient{k} = [1:m]; end
    end    
  end
  
  % Use comma-separated list syntax for N-D indexing.
  if useLinearShift,
    y = zeros(size(x));
    y(recipient{:}) = x(idx{:});
  else
    y = x(idx{:});
  end
return
%####################################################################
%####################################################################
%####################################################################
