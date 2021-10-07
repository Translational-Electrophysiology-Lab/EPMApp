function varargout = GUI_EGM_map_viewer_Ensite_modify_all_info(varargin)
% GUI_EGM_MAP_VIEWER_ENSITE_MODIFY_ALL_INFO MATLAB code for GUI_EGM_map_viewer_Ensite_modify_all_info.fig
%      GUI_EGM_MAP_VIEWER_ENSITE_MODIFY_ALL_INFO, by itself, creates a new GUI_EGM_MAP_VIEWER_ENSITE_MODIFY_ALL_INFO or raises the existing
%      singleton*.
%
%      H = GUI_EGM_MAP_VIEWER_ENSITE_MODIFY_ALL_INFO returns the handle to a new GUI_EGM_MAP_VIEWER_ENSITE_MODIFY_ALL_INFO or the handle to
%      the existing singleton*.
%
%      GUI_EGM_MAP_VIEWER_ENSITE_MODIFY_ALL_INFO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_EGM_MAP_VIEWER_ENSITE_MODIFY_ALL_INFO.M with the given input arguments.
%
%      GUI_EGM_MAP_VIEWER_ENSITE_MODIFY_ALL_INFO('Property','Value',...) creates a new GUI_EGM_MAP_VIEWER_ENSITE_MODIFY_ALL_INFO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_EGM_map_viewer_Ensite_modify_all_info_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_EGM_map_viewer_Ensite_modify_all_info_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_EGM_map_viewer_Ensite_modify_all_info

% Last Modified by GUIDE v2.5 19-May-2016 10:58:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_EGM_map_viewer_Ensite_modify_all_info_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_EGM_map_viewer_Ensite_modify_all_info_OutputFcn, ...
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


% --- Executes just before GUI_EGM_map_viewer_Ensite_modify_all_info is made visible.
function GUI_EGM_map_viewer_Ensite_modify_all_info_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_EGM_map_viewer_Ensite_modify_all_info (see VARARGIN)

% Choose default command line output for GUI_EGM_map_viewer_Ensite_modify_all_info
% handles.output = [];

handles.data = guidata(varargin{1});
%
dat = {0,120,100;180,350,25};
set(handles.uitable1,'data',dat,'ColumnEditable',[true(1,3)]);
handles.uitable1.RowName = {'AT','RT'};
handles.tag_do_output = 1;
ibs = handles.data.tag_beat_number.String;
set(handles.tag_edit_beats,'string',ibs)
set(handles.tag_button_beats_ALL,'value',0);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_EGM_map_viewer_Ensite_modify_all_info wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_EGM_map_viewer_Ensite_modify_all_info_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.uitable1.Data;
if ~handles.tag_do_output
    Output = [];
else
    Output.tab = handles.uitable1.Data;
    Output.do_AT = handles.tag_AT.Value;
    Output.do_RT = handles.tag_RT.Value;
    Output.do_type = handles.tag_type_of_mod.Value;
    str = handles.tag_edit_beats.String;
    i_u = find(str=='-');
    if isempty(i_u)
        Output.beats = str2double(str);
    else
        ib01 = str2double(str(1:i_u-1));
        ib02 = str2double(str(i_u+1:end));
        Output.beats = [ib01 : ib02];
    end
end
varargout{1} = Output;

delete(handles.figure1);



% --- Executes on button press in tag_AT.
function tag_AT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_AT


% --- Executes on button press in tag_RT.
function tag_RT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_RT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_RT


% --- Executes on selection change in tag_type_of_mod.
function tag_type_of_mod_Callback(hObject, eventdata, handles)
% hObject    handle to tag_type_of_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_type_of_mod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_type_of_mod


% --- Executes during object creation, after setting all properties.
function tag_type_of_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_type_of_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_OK.
function tag_OK_Callback(hObject, eventdata, handles)
% hObject    handle to tag_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mode_type = handles.tag_type_of_mod.Value;
do_AT = handles.tag_AT.Value;
do_RT = handles.tag_RT.Value;

AT_min = handles.uitable1.Data{1,1};
AT_max = handles.uitable1.Data{1,2};
RT_min = handles.uitable1.Data{2,1};
RT_max = handles.uitable1.Data{2,2};

Fcut_AT = handles.uitable1.Data{1,3};
Fcut_RT = handles.uitable1.Data{2,3};


delete(findall(gcf,'Type','hggroup'));
% handles.output = get(handles.uitable1,'data');
% guidata(hObject,handles)
% 
close(handles.figure1)

% --- Executes on button press in tag_cancel.
function tag_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to tag_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.tag_do_output = 0;
guidata(hObject,handles);
close(handles.figure1)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
% The GUI is still in UIWAIT, us UIRESUME
uiresume(hObject);
else
% The GUI is no longer waiting, just close it
delete(hObject);
end


% --- Executes on button press in tag_button_beats_ALL.
function tag_button_beats_ALL_Callback(hObject, eventdata, handles)
% hObject    handle to tag_button_beats_ALL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_button_beats_ALL
if handles.tag_button_beats_ALL.Value
    str = ['1-',num2str(length(handles.data.matfile.spikes))];
    set(handles.tag_edit_beats,'string',str)
else
    str = handles.data.tag_beat_number.String;
    set(handles.tag_edit_beats,'string',str)
end


function tag_edit_beats_Callback(hObject, eventdata, handles)
% hObject    handle to tag_edit_beats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_edit_beats as text
%        str2double(get(hObject,'String')) returns contents of tag_edit_beats as a double


% --- Executes during object creation, after setting all properties.
function tag_edit_beats_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_edit_beats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
