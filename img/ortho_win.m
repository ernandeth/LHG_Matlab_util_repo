function varargout = ortho_win(varargin)
% ORTHO_WIN M-file for ortho_win.fig
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
%      ORTHO_WIN, by itself, creates a new ORTHO_WIN or raises the existing
%      singleton*.
%
%      H = ORTHO_WIN returns the handle to a new ORTHO_WIN or the handle to
%      the existing singleton*.
%
%      ORTHO_WIN('Property','Value',...) creates a new ORTHO_WIN using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ortho_win_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      ORTHO_WIN('CALLBACK') and ORTHO_WIN('CALLBACK',hObject,...) call the
%      local function named CALLBACK in ORTHO_WIN.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Edit the above text to modify the response to help ortho_win

% Last Modified by GUIDE v2.5 02-Aug-2014 19:01:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ortho_win_OpeningFcn, ...
                   'gui_OutputFcn',  @ortho_win_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ortho_win is made visible.
function ortho_win_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ortho_win
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% if exist('~/matlab/img/LuisOrtho.wav')
%     snd = wavread('~/matlab/img/LuisOrtho.wav');
%     sound(snd,44100);
% end
% UIWAIT makes ortho_win wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% these are the boolean variables to determine if a field has been filled
% in or an option chosen ...

global args

    args.ROIsize = 0;
    args.ROItype = 'cube';
    args.threshold = 0;
    args.threshold2 = 0;
    args.onsets = [];
    args.window = 20;
    args.spm_file = [];
    args.spm_file2 = [];
    args.anat_file = [];
    args.tseries_path = [];
    args.tseries_file = [];
    args.tseries_file2 = [];
    args.doDetrend = 0;
    args.doGfilter = 0;
    args.doFFT = 0;    
    args.ignore_origin = 0;
    args.wscale = [];
    args.interact = 1;
    args.xyz=[];
    args.mask_file = [];
    args.output_name = 'Ortho';
    args.voxFile = [];
    args.doMovie = 0;
    args.causalMap = 0;
    args.doConnect = 0;
	args.doSufNec = 0;

% --- Outputs from this function are returned to the command line.
function varargout = ortho_win_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in GO_button.
function GO_button_Callback(hObject, eventdata, handles)
% hObject    handle to GO_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global args
    handle = findobj('Tag','AnatFile_field');
    args.anat_file = get(handle, 'String');
args

warning off

ortho2005(args)

return

% --- Executes on button press in STOP_button.
function STOP_button_Callback(hObject, eventdata, handles)
% hObject    handle to STOP_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global myfig EvFig histFig FFTfig ACTIVE

ACTIVE=0;

% Call local function for each fig type to find highest tag number
strMyFigLast = lochighestexistingfigure('\d{2}_myfig');
strEvFigLast = lochighestexistingfigure('\d{2}_evfig');
strHistFigLast = lochighestexistingfigure('\d{2}_histfig');
strFFTFigLast = lochighestexistingfigure('\d{2}_fftfig');

% Find highest tag number across all fig types and close all of that
% iteration.
casHighest = {strEvFigLast;strFFTFigLast;strHistFigLast;strMyFigLast};
casHighest = sort(casHighest);
strPatClose = casHighest{end};
hToClose = findobj('-regexp','tag',strPatClose);
close(hToClose);

return

% --- Executes during object creation, after setting all properties.
function ROItype_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROItype_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in ROItype_list.
function ROItype_list_Callback(hObject, eventdata, handles)
% hObject    handle to ROItype_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ROItype_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROItype_list
    global args
    ROItype = get(findobj('Tag','ROItype_list'), 'Value')
    switch ROItype
        case 1
            args.ROItype = 'cube'
        case 2
            args.ROItype = 'sphere'
            
        case 3
            args.ROItype = 'voxFile'
            
        case 4
            args.ROItype = 'maskFile'
    end
        

% --- Executes during object creation, after setting all properties.
function NN_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NN_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function NN_field_Callback(hObject, eventdata, handles)
% hObject    handle to NN_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NN_field as text
%        str2double(get(hObject,'String')) returns contents of NN_field as a double
global args
args.ROIsize = str2num(get(gco, 'String'))

% --- Executes on button press in AnatFile_button.
function AnatFile_button_Callback(hObject, eventdata, handles)
% hObject    handle to AnatFile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*.img;*.nii;*.nii.gz','Select  *.img;*.nii;*.nii.gz file');
    if isa(name,'double')
        return
    end
	name = strcat(path,name)
    handle = findobj('Tag','AnatFile_field');
    set(handle, 'String',name);
    args.anat_file = name;
return
    
    % --- Executes on button press in SM1_button.
function SM1_button_Callback(hObject, eventdata, handles)
% hObject    handle to SM1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*.img;*.nii;*.nii.gz','Select statistical *.img;*.nii;*.nii.gz file');
    if isa(name,'double')
        return
    end
    name = strcat(path,name)
    handle = findobj('Tag','SM1_field');
    set(handle, 'String',name);
    args.spm_file = name;
return
 

% --- Executes on button press in SM2_button.
function SM2_button_Callback(hObject, eventdata, handles)
% hObject    handle to SM2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*.img;*.nii;*.nii.gz','Select statistical *.img;*.nii;*.nii.gz file');
    if isa(name,'double')
        return
    end
	name = strcat(path,name)
    handle = findobj('Tag','SM2_field');
    set(handle, 'String',name);
    args.spm_file2 = name;
return

% --- Executes on button press in TS1_button.
function TS1_button_Callback(hObject, eventdata, handles)
% hObject    handle to TS1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*.img;*.nii;*.nii.gz','Select statistical *.img;*.nii;*.nii.gz file');
    if isa(name,'double')
        return
    end
	name = strcat(path,name)
    handle = findobj('Tag','TS1_field');
    set(handle, 'String',name);
    args.tseries_file = name;
return

% --- Executes on button press in TS2_button.
function TS2_button_Callback(hObject, eventdata, handles)
% hObject    handle to TS2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*.img;*.nii;*.nii.gz','Select statistical *.img;*.nii;*.nii.gz file');
	name = strcat(path,name)
    handle = findobj('Tag','TS2_field');
    set(handle, 'String',name);
    args.tseries2_file = name;
return


% --- Executes on button press in TS3_button.
function TS3_button_Callback(hObject, eventdata, handles)
% hObject    handle to TS3_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [name path] = uigetfile('*.img;*.nii;*.nii.gz','Select statistical *.img;*.nii;*.nii.gz file');
	name = strcat(path,name)
    handle = findobj('Tag','TS3_field');
    set(handle, 'String',name);

return


% --- Executes during object creation, after setting all properties.
function Th1_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Th1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Th1_field_Callback(hObject, eventdata, handles)
% hObject    handle to Th1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Th1_field as text
%        str2double(get(hObject,'String')) returns contents of Th1_field as a double
    global args
    str = get(findobj('Tag','Th1_field'), 'String');
    args.threshold = str2num(str);
    uiresume(gcbf);
% return
    
% --- Executes during object creation, after setting all properties.
function Th2_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Th2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Th2_field_Callback(hObject, eventdata, handles)
% hObject    handle to Th2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Th2_field as text
%        str2double(get(hObject,'String')) returns contents of Th2_field as a double
global args
    str = get(findobj('Tag','Th2_field'), 'String');
    args.threshold2 = str2num(str);
uiresume(gcbf)    
% return
    

% --- Executes during object creation, after setting all properties.
function ons1_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ons1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function ons1_field_Callback(hObject, eventdata, handles)
% hObject    handle to ons1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ons1_field as text
%        str2double(get(hObject,'String')) returns contents of ons1_field as a double
global args

str = get(hObject,'String')
if length(str)>0 
    args.onsets = str2num(str)
end
uiresume(gcbf);
% return

% --- Executes during object creation, after setting all properties.
function on2_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to on2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --- Executes during object creation, after setting all properties.
function ons3_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ons3_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ons3_field_Callback(hObject, eventdata, handles)
% hObject    handle to ons3_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ons3_field as text
%        str2double(get(hObject,'String')) returns contents of ons3_field as a double
    

% --- Executes during object creation, after setting all properties.
function TW_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TW_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function TW_field_Callback(hObject, eventdata, handles)
% hObject    handle to TW_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TW_field as text
%        str2double(get(hObject,'String')) returns contents of TW_field as a double
   global args
    str = get(hObject,'String')
    args.window = str2num(str);
%     return
uiresume(gcbf)

% --- Executes during object creation, after setting all properties.
function anat_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anat_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function SM2_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SM2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function SM2_field_Callback(hObject, eventdata, handles)
% hObject    handle to SM2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SM2_field as text
%        str2double(get(hObject,'String')) returns contents of SM2_field as a double
global args
args.spm_file2 = get(hObject,'String');

return

% --- Executes during object creation, after setting all properties.
function SM1_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SM1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function SM1_field_Callback(hObject, eventdata, handles)
% hObject    handle to SM1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SM1_field as text
%        str2double(get(hObject,'String')) returns contents of SM1_field as a double
global args
args.spm_file = get(hObject,'String');
return

% --- Executes during object creation, after setting all properties.
function TS3_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TS3_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function TS3_field_Callback(hObject, eventdata, handles)
% hObject    handle to TS3_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TS3_field as text
%        str2double(get(hObject,'String')) returns contents of TS3_field as a double



% --- Executes during object creation, after setting all properties.
function TS1_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TS1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function TS1_field_Callback(hObject, eventdata, handles)
% hObject    handle to TS1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TS1_field as text
%        str2double(get(hObject,'String')) returns contents of TS1_field as a double
   global args

    args.tseries_file = get(hObject,'String');

return

% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes on button press in VF_button.
function VF_button_Callback(hObject, eventdata, handles)
% hObject    handle to VF_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*.img;*.nii;*.nii.gz','Select  file with ROI');
	name = strcat(path,name)
    handle = findobj('Tag','VF_field');
    set(handle, 'String',name);
    if strcmp(args.ROItype,'maskFile')
        args.mask_file = name;
    else
        args.voxFile = name;
    end
    
return
    
% --- Executes during object creation, after setting all properties.
function VF_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VF_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function VF_field_Callback(hObject, eventdata, handles)
% hObject    handle to VF_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VF_field as text
%        str2double(get(hObject,'String')) returns contents of VF_field as a double
   global args
    str = get(hObject,'String')
    args.voxFile = str;

% --- Executes during object creation, after setting all properties.
function ons2_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ons2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ons2_field_Callback(hObject, eventdata, handles)
% hObject    handle to ons2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ons2_field as text
%        str2double(get(hObject,'String')) returns contents of ons2_field as a double

% --- Executes during object creation, after setting all properties.
function AnatFile_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnatFile_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function AnatFile_field_Callback(hObject, eventdata, handles)
% hObject    handle to AnatFile_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AnatFile_field as text
%        str2double(get(hObject,'String')) returns contents of AnatFile_field as a double
    global args
    str = get(hObject,'String')
    args.anat_file = str;

return


% --- Executes during object creation, after setting all properties.
function TS2_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TS2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function TS2_field_Callback(hObject, eventdata, handles)
% hObject    handle to TS2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TS2_field as text
%        str2double(get(hObject,'String')) returns contents of TS2_field as a double

return


% --- Executes on button press in detrend_cbx.
function detrend_cbx_Callback(hObject, eventdata, handles)
% hObject    handle to detrend_cbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detrend_cbx
global args
args.doDetrend = get(hObject,'Value');
return

% --- Executes on button press in GFilt_cbx.
function GFilt_cbx_Callback(hObject, eventdata, handles)
% hObject    handle to GFilt_cbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GFilt_cbx
global args
args.doGfilter = get(hObject,'Value');
return


% --- Executes on button press in FTbox.
function FTbox_Callback(hObject, eventdata, handles)
% hObject    handle to FTbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FTbox
global args
args.doFFT = get(hObject,'Value')
return



function output_name_Callback(hObject, eventdata, handles)
% hObject    handle to output_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_name as text
%        str2double(get(hObject,'String')) returns contents of output_name as a double
global args
args.output_name = get(hObject,'String');
return

% --- Executes during object creation, after setting all properties.
function output_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox10.
function Moviebox_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10
global args
args.doMovie = get(hObject,'Value')
return


% --- Executes on button press in checkbox11.
function CausalBox_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11
global args
args.causalMap = get(hObject,'Value')
return



function coordsBox_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double
global args

str = get(hObject, 'String')';
args.xyz = str2num(str');


% if ~isempty(str)
%     args.interact = 0;
% else
%     args.interact = 1;
% end
return

% --- Executes during object creation, after setting all properties.
function coordsBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox12.
function ConnectBox_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox12
global args
args.doConnect = get(hObject,'Value')
return



function wscale_field_Callback(hObject, eventdata, handles)
% hObject    handle to wscale_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wscale_field as text
%        str2double(get(hObject,'String')) returns contents of wscale_field as a double
global args
str = get(findobj('Tag','wscale_field'), 'String');
args.wscale = str2num(str)
uiresume(gcbf)
% return


% --- Executes during object creation, after setting all properties.
function wscale_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wscale_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over STOP_button.
function STOP_button_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to STOP_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on STOP_button and none of its controls.
function STOP_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to STOP_button (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function strHighest = lochighestexistingfigure(strPatTag)
% hExist = findobj('-regexp','tag','\d{2}_myfig');
hExist = findobj('-regexp','tag',strPatTag);
if isempty(hExist)
    strHighest = '';
else
    casTags = get(hExist,'tag');
    if ~ischar(casTags)
        casTags = sort(casTags);
        strTagHighest = casTags{end};
    else
        strTagHighest = casTags;
    end
    strHighest = strTagHighest(1:3);
end
        
% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13
global args
args.doSufNec = get(hObject,'Value')
return


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Th2_field.
function Th2_field_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Th2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
