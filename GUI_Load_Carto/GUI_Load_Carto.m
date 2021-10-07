function varargout = GUI_Load_Carto(varargin)
% GUI_LOAD_CARTO MATLAB code for GUI_Load_Carto.fig
%      GUI_LOAD_CARTO, by itself, creates a new GUI_LOAD_CARTO or raises the existing
%      singleton*.
%
%      H = GUI_LOAD_CARTO returns the handle to a new GUI_LOAD_CARTO or the handle to
%      the existing singleton*.
%
%      GUI_LOAD_CARTO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_LOAD_CARTO.M with the given input arguments.
%
%      GUI_LOAD_CARTO('Property','Value',...) creates a new GUI_LOAD_CARTO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Load_Carto_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Load_Carto_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Load_Carto

% Last Modified by GUIDE v2.5 21-May-2019 11:49:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_Load_Carto_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_Load_Carto_OutputFcn, ...
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


% --- Executes just before GUI_Load_Carto is made visible.
function GUI_Load_Carto_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Load_Carto (see VARARGIN)

% Choose default command line output for GUI_Load_Carto
handles.output = hObject;
pp = path;
if ~contains(pp,'GUI_egm_mFiles')
    d = pwd;ii=find(d=='\');addpath([d(1:ii(end)-1),'\GUI_egm_mFiles']);
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Load_Carto wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Load_Carto_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PATHNAME = uigetdir('','Select a directory'); % select a dir

Xcar = dir([PATHNAME,'\*_car.txt']);
car_names = cellfun(@(x) x(1:end-8),{Xcar.name},'un',0);
T = cell(length(car_names),6);
for i = 1:length(car_names)
    cc = LoadCartoCar_fun([PATHNAME,'\',Xcar(i).name],0);
    T{i,1} = car_names{i};
    if ~isempty(cc)
        % Get number of points
        Np = length(cc.Index_point);
        % Get type of catheter
        Xnavistar = dir([PATHNAME,'\',car_names{i},'*NAVISTAR_CONNECTOR_Eleclectrode_Positions_OnAnnotation.txt']);
        X20A = dir([PATHNAME,'\',car_names{i},'*MAGNETIC_20_POLE_A_CONNECTOR_Eleclectrode_Positions_OnAnnotation.txt']);
        X20B = dir([PATHNAME,'\',car_names{i},'*MAGNETIC_20_POLE_B_CONNECTOR_Eleclectrode_Positions_OnAnnotation.txt']);
        XCS = dir([PATHNAME,'\',car_names{i},'*CS_CONNECTOR_Eleclectrode_Positions_OnAnnotation.txt']);
        
        
        T{i,2} = Np;
        T{i,3} = length(Xnavistar);
        T{i,4} = length(X20A);
        T{i,5} = length(X20B);
        T{i,6} = length(XCS);
    else
        T(i,2:6) = num2cell(zeros(1,5));
    end
    
    
end
clear i

T = [num2cell([T{:,2}]>20)' T];

handles.uitable1.Data = T;
handles.mat.main_dir = PATHNAME;
guidata(handles.figure1,handles);

% --- Executes on button press in pushbutton_convert.
function pushbutton_convert_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_convert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tab = handles.uitable1.Data;
Tab(~[Tab{:,1}],:)=[];
do_pacing = handles.checkbox_pacing.Value;
do_markers = handles.checkbox_markers.Value;

ParamWin.do_window = handles.checkbox_Template.Value;
ParamWin.cc_QRS = handles.uitable_cc.Data{3,1};
ParamWin.cc_sig = handles.uitable_cc.Data{1,1};
ParamWin.cc_TW = handles.uitable_cc.Data{2,1};

for i_map = 1:size(Tab,1);
    filename = [handles.mat.main_dir,'\',Tab{i_map,2}];
    [signals,ParamSig] = Convert_Carto_fun(filename,do_pacing,ParamWin,do_markers);
    
end


% --- Executes on button press in checkbox_markers.
function checkbox_markers_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_markers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_markers


% --- Executes on button press in checkbox_pacing.
function checkbox_pacing_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_pacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_pacing


% --- Executes on button press in checkbox_Template.
function checkbox_Template_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Template
switch handles.checkbox_Template.Value
    case 0
        handles.uitable_cc.Data{1,1} = -1;
    case 1
        handles.uitable_cc.Data{1,1} = 0.8;
end


% --- Executes on button press in pushbutton_preview.
function pushbutton_preview_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dir_plot = handles.mat.main_dir;
GUI_plot_carto_point(dir_plot);