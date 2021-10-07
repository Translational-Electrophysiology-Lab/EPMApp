function varargout = ReLable_GUI(varargin)
% RELABLE_GUI MATLAB code for ReLable_GUI.fig
%      RELABLE_GUI, by itself, creates a new RELABLE_GUI or raises the existing
%      singleton*.
%
%      H = RELABLE_GUI returns the handle to a new RELABLE_GUI or the handle to
%      the existing singleton*.
%
%      RELABLE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RELABLE_GUI.M with the given input arguments.
%
%      RELABLE_GUI('Property','Value',...) creates a new RELABLE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ReLable_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ReLable_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ReLable_GUI

% Last Modified by GUIDE v2.5 21-Sep-2021 13:33:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ReLable_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ReLable_GUI_OutputFcn, ...
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


% --- Executes just before ReLable_GUI is made visible.
function ReLable_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ReLable_GUI (see VARARGIN)

% Choose default command line output for ReLable_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ReLable_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

handles.hObject_main = varargin{1};
guidata(hObject,handles)


% --- Outputs from this function are returned to the command line.
function varargout = ReLable_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
handles_main=guidata(handles.hObject_main);
dat =  cell(length(handles_main.ParamSig.Label),3);
dat(:,1)=handles_main.ParamSig.Label;
dat(:,2) = handles_main.ParamSig.Label;
dat(:,3) = {''};
set(handles.tag_table,'data',dat,'columneditable',[true true true],'columnformat',{'char','char','char'});




% --- Executes on button press in tag_done.
function tag_done_Callback(hObject, eventdata, handles)
% hObject    handle to tag_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles_main = guidata(handles.hObject_main);
dat = get(handles.tag_table,'data');
handles_main.ParamSig.Label = dat(:,2);
handles_main.ParamSig.Label_ori = dat(:,2);
handles_main.ParamSig.Notes = dat(:,3);
%% update
guidata(handles.hObject_main,handles_main);
close(handles.figure1)
