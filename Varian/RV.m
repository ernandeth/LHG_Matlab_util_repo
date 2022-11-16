function varargout = RV(varargin)

% RV M-file for RV.fig
%      RV, by itself, creates a new RV or raises the existing
%      singleton*.
%
%      H = RV returns the handle to a new RV or the handle to
%      the existing singleton*.
%
%      RV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RV.M with the given input arguments.
%
%      RV('Property','Value',...) creates a new RV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RV_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RV_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RV

% Last Modified by GUIDE v2.5 13-Mar-2010 11:21:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RV_OpeningFcn, ...
                   'gui_OutputFcn',  @RV_OutputFcn, ...
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


% --- Executes just before RV is made visible.
function RV_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RV (see VARARGIN)

% Choose default command line output for RV
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RV wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RV_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

fig = openfig(mfilename,'reuse');
   %a=imread('SimpsonGray.jpg');
   
   a=imread('VarianScreen.png');
  % b=imread('');
   %subplot(3,2,3)
   
   H_image1=imshow(a);




% --------------------------------------------------------------------
function Varian_Open_Callback(hObject, eventdata, handles)
% hObject    handle to Varian_Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_kspace_Callback(hObject, eventdata, handles)
% hObject    handle to Open_kspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global image_space final_kspace
global pathname final_kspace
nslice=1;
[final_kspace,np,ntraces,nblocks]= ReadVarStart(nslice);

% --------------------------------------------------------------------
function Open_images_Callback(hObject, eventdata, handles)
% hObject    handle to Open_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear image_space dim rank1 bits1 
clear global image_space final_kspace nblocks array ntraces seqcon final_kspace
global pathname
ReadVarFDF



% --------------------------------------------------------------------
function Open_Spec_Callback(hObject, eventdata, handles)
% hObject    handle to Open_Spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_spec_Callback(hObject, eventdata, handles)
% hObject    handle to Open_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global final_kspace image_space np ntraces nblocks
global pathname final_kspace
nspec=1
[final_kspace,np,ntraces,nblocks]= ReadVarStartSpec(nspec);





