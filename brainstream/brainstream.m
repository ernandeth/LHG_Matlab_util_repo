function varargout = brainstream(varargin)
% BRAINSTREAM M-file for brainstream.fig
%      BRAINSTREAM, by itself, creates a new BRAINSTREAM or raises the existing
%      singleton*.
%
%      H = BRAINSTREAM returns the handle to a new BRAINSTREAM or the handle to
%      the existing singleton*.
%
%      BRAINSTREAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRAINSTREAM.M with the given input arguments.
%
%      BRAINSTREAM('Property','Value',...) creates a new BRAINSTREAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fasl01_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to brainstream_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help brainstream

% Last Modified by GUIDE v2.5 29-Apr-2011 13:44:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @brainstream_OpeningFcn, ...
                   'gui_OutputFcn',  @brainstream_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

global args


% --- Executes just before brainstream is made visible.
function brainstream_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fasl01_working (see VARARGIN)

% Choose default command line output for fasl01_working
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fasl01_working wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = brainstream_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in go_button.
function go_button_Callback(hObject, eventdata, handles)
% hObject    handle to go_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global args
str = get(findobj('Tag','inputData_edit'),'String');
str2 = get(findobj('Tag','inputData_edit2'),'String');
str3 = get(findobj('Tag','inputData_edit3'),'String');
str4 = get(findobj('Tag','inputData_edit4'),'String');
str5 = get(findobj('Tag','inputData_edit5'),'String');
if ~isempty(str)
    
    args.inFile = str;  
    args.inFile2 = str2;
    args.inFile3 = str3;
    args.inFile4 = str4;
    args.inFile5 = str5;
    args.doDespike = get(findobj('Tag','despiker_cb'),'Value');
    args.doRecon =   get(findobj('Tag','SpiralRecon_cb'),'Value');

    args.doSliceTime = get(findobj('Tag','sliceTiming_cb'),'Value');
    args.doRealign = get(findobj('Tag','realignment_cb'),'Value');
    args.smoothSize= get(findobj('Tag','smoothing_cb'),'Value') *3 ;

    tmp=get(findobj('Tag','subtraction_cb'),'Value');
    if tmp
        args.isSubtracted=1;
        set(findobj('Tag','isSubtracted_cb'),'Value',1);
        args.subType = get(findobj('Tag','subtraction_lb'),'Value');
    else
        args.subType = 0;
    end

    args.physCorr = get(findobj('Tag','physio_cb'),'Value');
    args.physFile = get(findobj('Tag','physioName_edit'),'String');
    args.zflip = get(findobj('Tag', 'zflip'), 'Value');
    
    
    args.doGLM = get(findobj('Tag','doGLM_cb'),'Value');
    if args.doGLM ~= 0
        xname = get(findobj('Tag','loadDesignMatrix_edit'),'String');
        X = load(xname);
        c = eval(get(findobj('Tag','contrasts_edit'),'String'));
        
        args.contrasts = c;
        args.designMat = X;
        args.isSubtracted = get(findobj('Tag','isSubtracted_cb'),'Value');

        disp('Note Sizes of contrasts and Design Matrix:');
        whos c X
    end
    
    args.doQuant = get(findobj('Tag','doQuant_cb'),'Value');

    args.aslParms.TR = str2num(get(findobj('Tag','tr_edit'),'String'));
    args.aslParms.Ttag = str2num(get(findobj('Tag','Ttag_edit'),'String'));
    args.aslParms.Tdelay = str2num(get(findobj('Tag','Tdelay_edit'),'String'));
    args.aslParms.Ttransit = str2num(get(findobj('Tag','Ttransit_edit'),'String'));
    args.aslParms.inv_alpha = str2num(get(findobj('Tag','alpha_edit'),'String'));

    args.doLightbox = get(findobj('Tag','doLightbox_cb'),'Value');
    args.doAnantPreprocess = get(findobj('Tag', 'preprocess_anat'), 'Value');
    if args.doAnantPreprocess ~= 0
        anatdir = get(findobj('Tag', 'anatpathname'), 'String');
    end
    args.doOrtho = get(findobj('Tag','displayZmap_cb'),'Value');
    args.doCoreg = get(findobj('Tag','useOverlay_checkbox'),'Value');
    if args.doCoreg ~= 0
        overlayfile = get(findobj('Tag','overlayfile_edit'),'String');
        end
    args.useSPGR = get(findobj('Tag','useSPGR_checkbox'),'Value');
    if args.useSPGR ~= 0
        spgrfile = get(findobj('Tag','inputspgr_edit'),'String');
    end
    args.doTransfer = get(findobj('Tag', 'transfer'), 'Value');
    
    
    %%  Here is the call to the main function:  %%%%%%%
    asl_spm01(args);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    disp('NO DATA to WORK ON!')
end
% --- Executes on button press in inputData_button.
function inputData_button_Callback(hObject, eventdata, handles)
% hObject    handle to inputData_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str=uigetfile('*.*','Input Data (kpace OR image space)')
set(findobj('Tag','inputData_edit'),'String', str);

function inputData_edit_Callback(hObject, eventdata, handles)
% hObject    handle to inputData_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputData_edit as text
%        str2double(get(hObject,'String')) returns contents of inputData_edit as a double


% --- Executes during object creation, after setting all properties.
function inputData_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputData_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in inputData_button2.
function inputData_button2_Callback(hObject, eventdata, handles)
% hObject    handle to inputData_button2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str2=uigetfile('*.*','Input Data (kpace OR image space)')
set(findobj('Tag','inputData_edit2'),'String', str2);

function inputData_edit2_Callback(hObject, eventdata, handles)
% hObject    handle to inputData_edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputData_edit2 as text
%        str2double(get(hObject,'String')) returns contents of inputData_edit2 as a double


% --- Executes during object creation, after setting all properties.
function inputData_edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputData_edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in inputData_button3.
function inputData_button3_Callback(hObject, eventdata, handles)
% hObject    handle to inputData_button3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str3=uigetfile('*.*','Input Data (kpace OR image space)')
set(findobj('Tag','inputData_edit3'),'String', str3);

function inputData_edit3_Callback(hObject, eventdata, handles)
% hObject    handle to inputData_edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputData_edit3 as text
%        str2double(get(hObject,'String')) returns contents of inputData_edit3 as a double


% --- Executes during object creation, after setting all properties.
function inputData_edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputData_edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in inputData_button4.
function inputData_button4_Callback(hObject, eventdata, handles)
% hObject    handle to inputData_button4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str4=uigetfile('*.*','Input Data (kpace OR image space)')
set(findobj('Tag','inputData_edit4'),'String', str4);

function inputData_edit4_Callback(hObject, eventdata, handles)
% hObject    handle to inputData_edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputData_edit4 as text
%        str2double(get(hObject,'String')) returns contents of inputData_edit4 as a double


% --- Executes during object creation, after setting all properties.
function inputData_edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputData_edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in inputData_button5.
function inputData_button5_Callback(hObject, eventdata, handles)
% hObject    handle to inputData_button5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str5=uigetfile('*.*','Input Data (kpace OR image space)')
set(findobj('Tag','inputData_edit5'),'String', str5);


function inputData_edit5_Callback(hObject, eventdata, handles)
% hObject    handle to inputData_edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputData_edit5 as text
%        str2double(get(hObject,'String')) returns contents of inputData_edit5 as a double


% --- Executes during object creation, after setting all properties.
function inputData_edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputData_edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in SpiralRecon_cb.
function SpiralRecon_cb_Callback(hObject, eventdata, handles)
% hObject    handle to SpiralRecon_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SpiralRecon_cb
if get(gco,'Value')==1
   set(findobj('Tag', 'despiker_cb'), 'Enable','on');
else
   set(findobj('Tag', 'despiker_cb'), 'Enable','off');
end

% --- Executes on button press in despiker_cb.
function despiker_cb_Callback(hObject, eventdata, handles)
% hObject    handle to despiker_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of despiker_cb


% --- Executes on button press in realignment_cb.
function realignment_cb_Callback(hObject, eventdata, handles)
% hObject    handle to realignment_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of realignment_cb


% --- Executes on button press in smoothing_cb.
function smoothing_cb_Callback(hObject, eventdata, handles)
% hObject    handle to smoothing_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of smoothing_cb


% --- Executes on button press in subtraction_cb.
function subtraction_cb_Callback(hObject, eventdata, handles)
% hObject    handle to subtraction_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subtraction_cb
if get(gco,'Value')==1
   set(findobj('Tag', 'subtraction_lb'), 'Enable','on');
   set(findobj('Tag','isSubtracted_cb'),'Value',1);
else
   set(findobj('Tag', 'subtraction_lb'), 'Enable','off');
end

% --- Executes on button press in physio_cb.
function physio_cb_Callback(hObject, eventdata, handles)
% hObject    handle to physio_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of physio_cb
if get(gco,'Value')==1
   set(findobj('Tag', 'loadPhysioData_txt'), 'Enable','on');
   set(findobj('Tag', 'physioName_edit'), 'Enable','on');
else
   set(findobj('Tag', 'loadPhysioData_txt'), 'Enable','off');
   set(findobj('Tag', 'physioName_edit'), 'Enable','off');
   
end

% --- Executes on selection change in subtraction_lb.
function subtraction_lb_Callback(hObject, eventdata, handles)
% hObject    handle to subtraction_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns subtraction_lb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from subtraction_lb


% --- Executes during object creation, after setting all properties.
function subtraction_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subtraction_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadPhysioData_txt.
function loadPhysioData_txt_Callback(hObject, eventdata, handles)
% hObject    handle to loadPhysioData_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str=uigetfile('*.*','Load Physio Data File');
set(findobj('Tag', 'loadPhysioData_edit'), 'String',str);

function loadDesignMatrix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to loadDesignMatrix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadDesignMatrix_edit as text
%        str2double(get(hObject,'String')) returns contents of loadDesignMatrix_edit as a double


% --- Executes during object creation, after setting all properties.
function loadDesignMatrix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadDesignMatrix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


   
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in doGLM_cb.
function doGLM_cb_Callback(hObject, eventdata, handles)
% hObject    handle to doGLM_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doGLM_cb
if get(gco,'Value')==1
   set(findobj('Tag', 'buildDesignMatrix_button'), 'Enable','on');
   set(findobj('Tag', 'loadDesignMatrix_button'), 'Enable','on');
   set(findobj('Tag', 'loadDesignMatrix_edit'), 'Enable','on');
   set(findobj('Tag', 'contrasts_edit'), 'Enable','on');
   set(findobj('Tag', 'contrasts_text'), 'Enable','on');
   set(findobj('Tag', 'displayZmap_cb'), 'Enable','on');
   
else
   set(findobj('Tag', 'buildDesignMatrix_button'), 'Enable','off');
   set(findobj('Tag', 'loadDesignMatrix_button'), 'Enable','off');
   set(findobj('Tag', 'loadDesignMatrix_edit'), 'Enable','off');
   set(findobj('Tag', 'contrasts_edit'), 'Enable','off');
   set(findobj('Tag', 'contrasts_text'), 'Enable','off');
   set(findobj('Tag', 'displayZmap_cb'), 'Enable','off');
      
end

% --- Executes on button press in displayZmap_cb.
function displayZmap_cb_Callback(hObject, eventdata, handles)
% hObject    handle to displayZmap_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayZmap_cb


% --- Executes on button press in loadDesignMatrix_button.
function loadDesignMatrix_button_Callback(hObject, eventdata, handles)
% hObject    handle to loadDesignMatrix_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, pth]=uigetfile('*.*','Load Design Matrix File');
str = fullfile(pth,file);
set(findobj('Tag', 'loadDesignMatrix_edit'), 'String',str);
X=load(str);
figure; imagesc(X); colormap gray ; title('Design Matrix');

function physioName_edit_Callback(hObject, eventdata, handles)
% hObject    handle to physioName_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of physioName_edit as text
%        str2double(get(hObject,'String')) returns contents of physioName_edit as a double


% --- Executes during object creation, after setting all properties.
function physioName_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to physioName_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in doQuant_cb.
function doQuant_cb_Callback(hObject, eventdata, handles)
% hObject    handle to doQuant_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doQuant_cb
if get(gco,'Value')==1
   set(findobj('Tag', 'doLightbox_cb'), 'Enable','on');
else
   set(findobj('Tag', 'doLightbox_cb'), 'Enable','off');
end

% --- Executes on button press in buildDesignMatrix_button.
function buildDesignMatrix_button_Callback(hObject, eventdata, handles)
% hObject    handle to buildDesignMatrix_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ASLmatBuilder01

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function contrasts_edit_Callback(hObject, eventdata, handles)
% hObject    handle to contrasts_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contrasts_edit as text
%        str2double(get(hObject,'String')) returns contents of contrasts_edit as a double


% --- Executes during object creation, after setting all properties.
function contrasts_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrasts_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_edit_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_edit as text
%        str2double(get(hObject,'String')) returns contents of alpha_edit as a double


% --- Executes during object creation, after setting all properties.
function alpha_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ttransit_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Ttransit_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ttransit_edit as text
%        str2double(get(hObject,'String')) returns contents of Ttransit_edit as a double


% --- Executes during object creation, after setting all properties.
function Ttransit_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ttransit_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tr_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tr_edit as text
%        str2double(get(hObject,'String')) returns contents of tr_edit as a double


% --- Executes during object creation, after setting all properties.
function tr_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TR_edit_Callback(hObject, eventdata, handles)
% hObject    handle to TR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TR_edit as text
%        str2double(get(hObject,'String')) returns contents of TR_edit as a double


% --- Executes during object creation, after setting all properties.
function TR_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ttag_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Ttag_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ttag_edit as text
%        str2double(get(hObject,'String')) returns contents of Ttag_edit as a double


% --- Executes during object creation, after setting all properties.
function Ttag_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ttag_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tdelay_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Tdelay_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tdelay_edit as text
%        str2double(get(hObject,'String')) returns contents of Tdelay_edit as a double


% --- Executes during object creation, after setting all properties.
function Tdelay_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tdelay_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sliceTiming_cb.
function sliceTiming_cb_Callback(hObject, eventdata, handles)
% hObject    handle to sliceTiming_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sliceTiming_cb


% --- Executes on button press in doLightbox_cb.
function doLightbox_cb_Callback(hObject, eventdata, handles)
% hObject    handle to doLightbox_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doLightbox_cb


% --- Executes on button press in isSubtracted_cb.
function isSubtracted_cb_Callback(hObject, eventdata, handles)
% hObject    handle to isSubtracted_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isSubtracted_cb


% --- Executes on button press in useOverlay_checkbox.
function useOverlay_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to useOverlay_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useOverlay_checkbox


% --- Executes on button press in InputOverlayFile.
function InputOverlayFile_Callback(hObject, eventdata, handles)
% hObject    handle to InputOverlayFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, pth]=uigetfile('*.*','Input Overlay File');
str = fullfile(pth,file);
set(findobj('Tag', 'overlayfile_edit'), 'String',str);


% --- Executes on button press in InputSPGRFile.
function InputSPGRFile_Callback(hObject, eventdata, handles)
% hObject    handle to InputSPGRFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, pth]=uigetfile('*.*','Input SPGR File');
str = fullfile(pth,file);
set(findobj('Tag', 'inputspgr_edit'), 'String',str);


function overlayfile_edit_Callback(hObject, eventdata, handles)
% hObject    handle to overlayfile_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of overlayfile_edit as text
%        str2double(get(hObject,'String')) returns contents of overlayfile_edit as a double


% --- Executes during object creation, after setting all properties.
function overlayfile_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overlayfile_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inputspgr_edit_Callback(hObject, eventdata, handles)
% hObject    handle to inputspgr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputspgr_edit as text
%        str2double(get(hObject,'String')) returns contents of inputspgr_edit as a double


% --- Executes during object creation, after setting all properties.
function inputspgr_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputspgr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in useSPGR_checkbox.
function useSPGR_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to useSPGR_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useSPGR_checkbox


% --- Executes on button press in zflip.
function zflip_Callback(hObject, eventdata, handles)
% hObject    handle to zflip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zflip


% --- Executes on button press in transfer.
function transfer_Callback(hObject, eventdata, handles)
% hObject    handle to transfer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of transfer


% --- Executes on button press in doAnatPreprocess.
function doAnatPreprocess_Callback(hObject, eventdata, handles)
% hObject    handle to doAnatPreprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doAnatPreprocess


% --- Executes on button press in anatpath.
function anatpath_Callback(hObject, eventdata, handles)
% hObject    handle to anatpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pth = uigetdir('*/*','Input path to anatomical folder');
str = path(pth);
set(findobj('Tag', 'anatpathname'), 'String',str);



function anatpathname_Callback(hObject, eventdata, handles)
% hObject    handle to anatpathname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of anatpathname as text
%        str2double(get(hObject,'String')) returns contents of anatpathname as a double

% --- Executes during object creation, after setting all properties.
function anatpathname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anatpathname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in input_file2transfer1.
function input_file2transfer1_Callback(hObject, eventdata, handles)
% hObject    handle to input_file2transfer1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function transfer_structural_Callback(hObject, eventdata, handles)
% hObject    handle to transfer_structural (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transfer_structural as text
%        str2double(get(hObject,'String')) returns contents of transfer_structural as a double


% --- Executes during object creation, after setting all properties.
function transfer_structural_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transfer_structural (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in file2transfer2.
function file2transfer2_Callback(hObject, eventdata, handles)
% hObject    handle to file2transfer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function transfer_contrast_Callback(hObject, eventdata, handles)
% hObject    handle to transfer_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transfer_contrast as text
%        str2double(get(hObject,'String')) returns contents of transfer_contrast as a double


% --- Executes during object creation, after setting all properties.
function transfer_contrast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transfer_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
