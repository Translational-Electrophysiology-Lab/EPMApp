function varargout = Alternans_Variability_Info(varargin)
% ALTERNANS_VARIABILITY_INFO MATLAB code for Alternans_Variability_Info.fig
%      ALTERNANS_VARIABILITY_INFO, by itself, creates a new ALTERNANS_VARIABILITY_INFO or raises the existing
%      singleton*.
%
%      H = ALTERNANS_VARIABILITY_INFO returns the handle to a new ALTERNANS_VARIABILITY_INFO or the handle to
%      the existing singleton*.
%
%      ALTERNANS_VARIABILITY_INFO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALTERNANS_VARIABILITY_INFO.M with the given input arguments.
%
%      ALTERNANS_VARIABILITY_INFO('Property','Value',...) creates a new ALTERNANS_VARIABILITY_INFO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Alternans_Variability_Info_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Alternans_Variability_Info_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Alternans_Variability_Info

% Last Modified by GUIDE v2.5 09-Jun-2014 14:50:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Alternans_Variability_Info_OpeningFcn, ...
                   'gui_OutputFcn',  @Alternans_Variability_Info_OutputFcn, ...
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


% --- Executes just before Alternans_Variability_Info is made visible.
function Alternans_Variability_Info_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Alternans_Variability_Info (see VARARGIN)

% Choose default command line output for Alternans_Variability_Info
handles.output = hObject;
handles.hObject_main = varargin{1};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Alternans_Variability_Info wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Alternans_Variability_Info_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = [];

% --- Executes on button press in tag_do_ARIA.
function tag_do_ARIA_Callback(hObject, eventdata, handles)
% hObject    handle to tag_do_ARIA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_do_ARIA


% --- Executes on button press in tag_do_RTA.
function tag_do_RTA_Callback(hObject, eventdata, handles)
% hObject    handle to tag_do_RTA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_do_RTA



function tag_n_surr_GLRT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_n_surr_GLRT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_n_surr_GLRT as text
%        str2double(get(hObject,'String')) returns contents of tag_n_surr_GLRT as a double


% --- Executes during object creation, after setting all properties.
function tag_n_surr_GLRT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_n_surr_GLRT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_n_surr_SM_Callback(hObject, eventdata, handles)
% hObject    handle to tag_n_surr_SM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_n_surr_SM as text
%        str2double(get(hObject,'String')) returns contents of tag_n_surr_SM as a double


% --- Executes during object creation, after setting all properties.
function tag_n_surr_SM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_n_surr_SM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_do_DTA.
function tag_do_DTA_Callback(hObject, eventdata, handles)
% hObject    handle to tag_do_DTA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_do_DTA


% --- Executes on button press in tag_do_TpA.
function tag_do_TpA_Callback(hObject, eventdata, handles)
% hObject    handle to tag_do_TpA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_do_TpA


% --- Executes on button press in tag_ok.
function tag_ok_Callback(hObject, eventdata, handles)
% hObject    handle to tag_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Alt_input.do_ARIA = get(handles.tag_do_ARIA,'value');
Alt_input.do_RTA = get(handles.tag_do_RTA,'value');
Alt_input.do_DTA = get(handles.tag_do_DTA,'value');
Alt_input.do_TpA = get(handles.tag_do_TpA,'value');

Alt_input.GLRT_n_surr = str2double(get(handles.tag_n_surr_GLRT,'string'));
Alt_input.SM_n_surr = str2double(get(handles.tag_n_surr_SM,'string'));

Alt_input.do_entire_series = get(handles.tag_entire_series,'value');
Alt_input.do_divide_by_CL = get(handles.tag_divide_by_CL,'value');

%%
data = guidata(handles.hObject_main);
data.Alt_input = Alt_input;
% update
guidata(handles.hObject_main,data);
close(handles.figure1)

% --- Executes on button press in tag_cancel.
function tag_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to tag_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = guidata(handles.hObject_main);
data.Alt_input = [];
% update
guidata(handles.hObject_main,data);
close(handles.figure1)


% --- Executes on button press in tag_entire_series.
function tag_entire_series_Callback(hObject, eventdata, handles)
% hObject    handle to tag_entire_series (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_entire_series


% --- Executes on button press in tag_divide_by_CL.
function tag_divide_by_CL_Callback(hObject, eventdata, handles)
% hObject    handle to tag_divide_by_CL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_divide_by_CL
