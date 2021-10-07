function varargout = Markers_input_ECGI(varargin)
% MARKERS_INPUT_ECGI MATLAB code for Markers_input_ECGI.fig
%      MARKERS_INPUT_ECGI, by itself, creates a new MARKERS_INPUT_ECGI or raises the existing
%      singleton*.
%
%      H = MARKERS_INPUT_ECGI returns the handle to a new MARKERS_INPUT_ECGI or the handle to
%      the existing singleton*.
%
%      MARKERS_INPUT_ECGI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARKERS_INPUT_ECGI.M with the given input arguments.
%
%      MARKERS_INPUT_ECGI('Property','Value',...) creates a new MARKERS_INPUT_ECGI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Markers_input_ECGI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Markers_input_ECGI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Markers_input_ECGI

% Last Modified by GUIDE v2.5 01-Feb-2016 13:49:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Markers_input_ECGI_OpeningFcn, ...
    'gui_OutputFcn',  @Markers_input_ECGI_OutputFcn, ...
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


% --- Executes just before Markers_input_ECGI is made visible.
function Markers_input_ECGI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Markers_input_ECGI (see VARARGIN)

data = guidata(varargin{1});
handles.hObject_main = data.figure1;
% Choose default command line output for Markers_input_ECGI
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Markers_input_ECGI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Markers_input_ECGI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = [];








function tag_AT_text_Callback(hObject, eventdata, handles)
% hObject    handle to tag_AT_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_AT_text as text
%        str2double(get(hObject,'String')) returns contents of tag_AT_text as a double


% --- Executes during object creation, after setting all properties.
function tag_AT_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_AT_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_filter_text_Callback(hObject, eventdata, handles)
% hObject    handle to tag_filter_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_filter_text as text
%        str2double(get(hObject,'String')) returns contents of tag_filter_text as a double


% --- Executes during object creation, after setting all properties.
function tag_filter_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_filter_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function tag_RT_HFcut_Callback(hObject, eventdata, handles)
% hObject    handle to tag_RT_HFcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_RT_HFcut as text
%        str2double(get(hObject,'String')) returns contents of tag_RT_HFcut as a double


% --- Executes during object creation, after setting all properties.
function tag_RT_HFcut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_RT_HFcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_RT_LFcut_Callback(hObject, eventdata, handles)
% hObject    handle to tag_RT_HFcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_RT_HFcut as text
%        str2double(get(hObject,'String')) returns contents of tag_RT_HFcut as a double


% --- Executes during object creation, after setting all properties.
function tag_RT_LFcut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_RT_HFcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_ATmax_Callback(hObject, eventdata, handles)
% hObject    handle to tag_ATmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_ATmax as text
%        str2double(get(hObject,'String')) returns contents of tag_ATmax as a double


% --- Executes during object creation, after setting all properties.
function tag_ATmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_ATmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_RTmax_Callback(hObject, eventdata, handles)
% hObject    handle to tag_RTmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_RTmax as text
%        str2double(get(hObject,'String')) returns contents of tag_RTmax as a double


% --- Executes during object creation, after setting all properties.
function tag_RTmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_RTmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_RTmin_Callback(hObject, eventdata, handles)
% hObject    handle to tag_RTmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_RTmin as text
%        str2double(get(hObject,'String')) returns contents of tag_RTmin as a double


% --- Executes during object creation, after setting all properties.
function tag_RTmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_RTmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in tag_OK.
function tag_OK_Callback(hObject, eventdata, handles)
% hObject    handle to tag_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MarkersInput.RT_BW = [str2double(get(handles.tag_RT_HFcut,'string')) str2double(get(handles.tag_RT_LFcut,'string'))];
MarkersInput.AT_BW = [str2double(get(handles.tag_AT_HFcut,'string')) str2double(get(handles.tag_AT_LFcut,'string'))];

MarkersInput.DTmax = str2double(get(handles.tag_ATmax,'string'));
MarkersInput.min_RT = str2double(get(handles.tag_RTmin,'string'));
MarkersInput.max_RT = str2double(get(handles.tag_RTmax,'string'));

%%
data = guidata(handles.hObject_main);
data.MarkersInput = MarkersInput;
% update
guidata(handles.hObject_main,data);


close(handles.figure1)




function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_AT_LFcut_Callback(hObject, eventdata, handles)
% hObject    handle to tag_AT_LFcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_AT_LFcut as text
%        str2double(get(hObject,'String')) returns contents of tag_AT_LFcut as a double


% --- Executes during object creation, after setting all properties.
function tag_AT_LFcut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_AT_LFcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_AT_HFcut_Callback(hObject, eventdata, handles)
% hObject    handle to tag_AT_HFcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_AT_HFcut as text
%        str2double(get(hObject,'String')) returns contents of tag_AT_HFcut as a double


% --- Executes during object creation, after setting all properties.
function tag_AT_HFcut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_AT_HFcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_RT_text_Callback(hObject, eventdata, handles)
% hObject    handle to tag_RT_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_RT_text as text
%        str2double(get(hObject,'String')) returns contents of tag_RT_text as a double


% --- Executes during object creation, after setting all properties.
function tag_RT_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_RT_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
