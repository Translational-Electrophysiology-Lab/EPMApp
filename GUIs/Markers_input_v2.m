function varargout = Markers_input_v2(varargin)
% MARKERS_INPUT_V2 MATLAB code for Markers_input_v2.fig
%      MARKERS_INPUT_V2, by itself, creates a new MARKERS_INPUT_V2 or raises the existing
%      singleton*.
%
%      H = MARKERS_INPUT_V2 returns the handle to a new MARKERS_INPUT_V2 or the handle to
%      the existing singleton*.
%
%      MARKERS_INPUT_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARKERS_INPUT_V2.M with the given input arguments.
%
%      MARKERS_INPUT_V2('Property','Value',...) creates a new MARKERS_INPUT_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Markers_input_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Markers_input_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Markers_input_v2

% Last Modified by GUIDE v2.5 14-Jul-2018 16:08:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Markers_input_v2_OpeningFcn, ...
    'gui_OutputFcn',  @Markers_input_v2_OutputFcn, ...
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


% --- Executes just before Markers_input_v2 is made visible.
function Markers_input_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Markers_input_v2 (see VARARGIN)

% Choose default command line output for Markers_input_v2
handles.output = hObject;
handles.hObject_main = varargin{1};

handles_main = guidata(handles.hObject_main);
handles_main.MarkersInput.ok = 0;
guidata(handles.hObject_main,handles_main);

% if handles_main.ParamSig.frequency>5000
if isfield(handles_main.ParamSig,'geoname')
    if isequal(handles_main.ParamSig.geoname,'MEA');
        set(handles.tag_LFcut,'string','5')
        set(handles.tag_HFcut,'string','400')
        set(handles.tag_ATmax,'string','15')
        set(handles.tag_RTmin,'string','20')
        set(handles.tag_RTmax,'string','80')
        
        
        set(handles.tag_SB_min,'string','5')
        set(handles.tag_SB_max,'string','400')
        set(handles.tag_NB_min,'string','400')
        set(handles.tag_NB_max,'string','1000')
        
        set(handles.tag_S1_resti,'string',' ')
        set(handles.tag_width_spike_replace,'string','1')
    end
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Markers_input_v2 wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Markers_input_v2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = [];


% --- Executes on selection change in tag_Type_DTRT.
function tag_Type_DTRT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Type_DTRT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_Type_DTRT contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_Type_DTRT
if handles.tag_Type_DTRT.Value==2
    set(handles.tag_LFcut,'string','0.5')
    set(handles.tag_HFcut,'string','25')
    set(handles.tag_LFcut_AT,'string','0')
    set(handles.tag_HFcut_AT,'string','100')
    
    set(handles.tag_ATmax,'string','100')
    set(handles.tag_RTmin,'string','150')
    set(handles.tag_RTmax,'string','350')
    set(handles.tag_DeltaT,'string','40')
    set(handles.tag_SB_min,'string','2')
    set(handles.tag_SB_max,'string','40')
    set(handles.tag_NB_min,'string','40')
    set(handles.tag_NB_max,'string','100')
    set(handles.tag_Thresh_spikes,'string',' ')
    set(handles.tag_width_spike_replace,'string',' ')
    
elseif handles.tag_Type_DTRT.Value==1
    set(handles.tag_LFcut,'string','0.5')
    set(handles.tag_HFcut,'string','25')
    set(handles.tag_LFcut_AT,'string','0')
    set(handles.tag_HFcut_AT,'string','100')
    set(handles.tag_DeltaT,'string','5')
    
    set(handles.tag_ATmax,'string','150')
    set(handles.tag_RTmin,'string','200')
    set(handles.tag_RTmax,'string','480')
    set(handles.tag_SB_min,'string','2')
    set(handles.tag_SB_max,'string','40')
    set(handles.tag_NB_min,'string','40')
    set(handles.tag_NB_max,'string','100')
    set(handles.tag_width_spike_replace,'string','15')
    
end


% --- Executes during object creation, after setting all properties.
function tag_Type_DTRT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_Type_DTRT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function tag_S_B_text_Callback(hObject, eventdata, handles)
% hObject    handle to tag_S_B_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_S_B_text as text
%        str2double(get(hObject,'String')) returns contents of tag_S_B_text as a double


% --- Executes during object creation, after setting all properties.
function tag_S_B_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_S_B_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_ARI_text_Callback(hObject, eventdata, handles)
% hObject    handle to tag_ARI_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_ARI_text as text
%        str2double(get(hObject,'String')) returns contents of tag_ARI_text as a double


% --- Executes during object creation, after setting all properties.
function tag_ARI_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_ARI_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_noise_B_text_Callback(hObject, eventdata, handles)
% hObject    handle to tag_noise_B_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_noise_B_text as text
%        str2double(get(hObject,'String')) returns contents of tag_noise_B_text as a double


% --- Executes during object creation, after setting all properties.
function tag_noise_B_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_noise_B_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_HFcut_Callback(hObject, eventdata, handles)
% hObject    handle to tag_HFcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_HFcut as text
%        str2double(get(hObject,'String')) returns contents of tag_HFcut as a double


% --- Executes during object creation, after setting all properties.
function tag_HFcut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_HFcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_LFcut_Callback(hObject, eventdata, handles)
% hObject    handle to tag_LFcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_LFcut as text
%        str2double(get(hObject,'String')) returns contents of tag_LFcut as a double


% --- Executes during object creation, after setting all properties.
function tag_LFcut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_LFcut (see GCBO)
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



function tag_SB_max_Callback(hObject, eventdata, handles)
% hObject    handle to tag_SB_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_SB_max as text
%        str2double(get(hObject,'String')) returns contents of tag_SB_max as a double


% --- Executes during object creation, after setting all properties.
function tag_SB_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_SB_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_SB_min_Callback(hObject, eventdata, handles)
% hObject    handle to tag_SB_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_SB_min as text
%        str2double(get(hObject,'String')) returns contents of tag_SB_min as a double


% --- Executes during object creation, after setting all properties.
function tag_SB_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_SB_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_NB_max_Callback(hObject, eventdata, handles)
% hObject    handle to tag_NB_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_NB_max as text
%        str2double(get(hObject,'String')) returns contents of tag_NB_max as a double


% --- Executes during object creation, after setting all properties.
function tag_NB_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_NB_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_NB_min_Callback(hObject, eventdata, handles)
% hObject    handle to tag_NB_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_NB_min as text
%        str2double(get(hObject,'String')) returns contents of tag_NB_min as a double


% --- Executes during object creation, after setting all properties.
function tag_NB_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_NB_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_DTRT_only.
function tag_DTRT_only_Callback(hObject, eventdata, handles)
% hObject    handle to tag_DTRT_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_DTRT_only


% --- Executes on button press in tag_OK.
function tag_OK_Callback(hObject, eventdata, handles)
% hObject    handle to tag_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = guidata(handles.hObject_main);
MarkersInput = data.MarkersInput;
MarkersInput.ok = 1;
MarkersInput.BW = [str2double(get(handles.tag_LFcut,'string')) str2double(get(handles.tag_HFcut,'string'))];
MarkersInput.BW_AT = [str2double(get(handles.tag_LFcut_AT,'string')) str2double(get(handles.tag_HFcut_AT,'string'))];
MarkersInput.DTmax = str2double(get(handles.tag_ATmax,'string'));
MarkersInput.min_RT = str2double(get(handles.tag_RTmin,'string'));
MarkersInput.max_RT = str2double(get(handles.tag_RTmax,'string'));
MarkersInput.do_DTRT_only_fast = get(handles.tag_DTRT_only,'value');
MarkersInput.Twidth_replace_spk = str2double(get(handles.tag_width_spike_replace,'string'));
MarkersInput.Thresh_replace_spk = str2double(get(handles.tag_Thresh_spikes,'string'));
MarkersInput.do_notch = handles.tag_do_notch.Value;
MarkersInput.DeltaT = str2double(get(handles.tag_DeltaT,'string'));

handles.SNR_Bsig = [str2double(get(handles.tag_SB_min,'string')),str2double(get(handles.tag_SB_max,'string'))];
handles.SNR_Bnoise = [str2double(get(handles.tag_NB_min,'string')),str2double(get(handles.tag_NB_max,'string'))];
str = get(handles.tag_Type_DTRT,'string');
handles.type_analysis = str{get(handles.tag_Type_DTRT,'value')};
%%
data.MarkersInput = MarkersInput;
data.type_analysis = handles.type_analysis;
data.SNR_Bsig = handles.SNR_Bsig;
data.SNR_Bnoise = handles.SNR_Bnoise;
% update
handles.hObject_main
guidata(handles.hObject_main,data);
close(handles.figure1)




function txt_spk_Callback(hObject, eventdata, handles)
% hObject    handle to txt_spk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_spk as text
%        str2double(get(hObject,'String')) returns contents of txt_spk as a double


% --- Executes during object creation, after setting all properties.
function txt_spk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_spk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_width_spike_replace_Callback(hObject, eventdata, handles)
% hObject    handle to tag_width_spike_replace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_width_spike_replace as text
%        str2double(get(hObject,'String')) returns contents of tag_width_spike_replace as a double


% --- Executes during object creation, after setting all properties.
function tag_width_spike_replace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_width_spike_replace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_Thresh_spikes_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Thresh_spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_Thresh_spikes as text
%        str2double(get(hObject,'String')) returns contents of tag_Thresh_spikes as a double


% --- Executes during object creation, after setting all properties.
function tag_Thresh_spikes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_Thresh_spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function tag_HFcut_AT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_HFcut_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_HFcut_AT as text
%        str2double(get(hObject,'String')) returns contents of tag_HFcut_AT as a double


% --- Executes during object creation, after setting all properties.
function tag_HFcut_AT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_HFcut_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_LFcut_AT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_LFcut_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_LFcut_AT as text
%        str2double(get(hObject,'String')) returns contents of tag_LFcut_AT as a double


% --- Executes during object creation, after setting all properties.
function tag_LFcut_AT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_LFcut_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_do_notch.
function tag_do_notch_Callback(hObject, eventdata, handles)
% hObject    handle to tag_do_notch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_do_notch



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



function tag_DeltaT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_DeltaT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_DeltaT as text
%        str2double(get(hObject,'String')) returns contents of tag_DeltaT as a double


% --- Executes during object creation, after setting all properties.
function tag_DeltaT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_DeltaT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_select.
function tag_select_Callback(hObject, eventdata, handles)
% hObject    handle to tag_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = guidata(handles.hObject_main);
maxCL = round(700/1000*data.ParamSig.frequency);
spikes = round((data.spikes - str2double(handles.tag_DeltaT.String))/1000*data.ParamSig.frequency);
L=min([max(diff(spikes)),maxCL]);
X = nan(L,length(spikes),size(data.S,2));

if length(spikes)>1
    for i= 1 : length(spikes)-1
        H = spikes(i) : spikes(i+1)-1;
        H(L:end)=[];H(H>size(data.S,1))=[];
        if H(1)<1
            
            X(-min(H)+H(H>0),i,:) =  data.S(H(H>0),:);
        else
            H(H<1)=[];
            X(1:length(H),i,:) = data.S(H,:);
        end
    end
    H = spikes(i+1):size(data.S,1);
    H(L:end)=[];
    X(1:length(H),i+1,:) = data.S(H,:);
else
    i = 0;
    H = spikes(i+1):size(data.S,1);
    H(L:end)=[];H(H<1)=[];
    X(1:length(H),i+1,:) = data.S(H,:);
    
end

GUI_Select_Windows_egm(X);
