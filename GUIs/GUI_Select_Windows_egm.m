function varargout = GUI_Select_Windows_egm(varargin)
% GUI_SELECT_WINDOWS_EGM MATLAB code for GUI_Select_Windows_egm.fig
%      GUI_SELECT_WINDOWS_EGM, by itself, creates a new GUI_SELECT_WINDOWS_EGM or raises the existing
%      singleton*.
%
%      H = GUI_SELECT_WINDOWS_EGM returns the handle to a new GUI_SELECT_WINDOWS_EGM or the handle to
%      the existing singleton*.
%
%      GUI_SELECT_WINDOWS_EGM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SELECT_WINDOWS_EGM.M with the given input arguments.
%
%      GUI_SELECT_WINDOWS_EGM('Property','Value',...) creates a new GUI_SELECT_WINDOWS_EGM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Select_Windows_egm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Select_Windows_egm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Select_Windows_egm

% Last Modified by GUIDE v2.5 08-Nov-2017 17:39:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_Select_Windows_egm_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_Select_Windows_egm_OutputFcn, ...
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


% --- Executes just before GUI_Select_Windows_egm is made visible.
function GUI_Select_Windows_egm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Select_Windows_egm (see VARARGIN)

% Choose default command line output for GUI_Select_Windows_egm
handles.output = hObject;
handles.S = varargin{1};
handles.Windows_template = cell(2,2);

handles.tag_slider.SliderStep = [1 10]/size(handles.S,2);
handles.tag_slider.Max = size(handles.S,2);
handles.tag_slider.Min = 1;
handles.tag_slider.Value = 1;

handles.listbox1.String = [repmat('chan-',[size(handles.S,3),1]) num2str([1:size(handles.S,3)]')];
handles.listbox1.Max = size(handles.S,3);

ff = findobj(0,'name','Markers_GUI_v3');
data = guidata(ff(1));
handles.ParamSig = data.ParamSig;
handles.spikes = data.spikes;

plot(handles.axes2,data.spikes(2:end)/1000,diff(handles.spikes),'.-');
if length(data.spikes)>1
    set(handles.axes2,'xlim',[0 data.spikes(end)]/1000);
end
xlabel(handles.axes2,'Time (s)')
tab_plot_button_Callback(hObject,[], handles);
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes GUI_Select_Windows_egm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Select_Windows_egm_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in tag_select_button.
function tag_select_button_Callback(hObject, eventdata, handles)
% hObject    handle to tag_select_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ax = [handles.axes1];
ibeat = str2double(handles.tag_point_number.String);
xx = get(ax(1),'xlim');
hold(ax(1),'on')
h(1) = patch([xx(1) xx(2) xx(2) xx(1)],[handles.axes1.YLim(1) handles.axes1.YLim(1) handles.axes1.YLim(2) handles.axes1.YLim(2)],[0 .8 .8],'parent',ax(1));
set(h,'FaceAlpha',0.1);
set(h,'EdgeColor','k');
set(h,'LineWidth',2);

hold(handles.axes2,'on');
rr = diff(handles.spikes);
if isempty(rr);rr=nan;end
a= plot(handles.axes2,handles.spikes(ibeat)/1000,rr(ibeat),'or','markerfacecolor','k');

%
switch handles.tag_select_type.Value;
    case 1
        set(h,'facecolor',[0 .8 .8],'EdgeColor',[0 .8 .8]);
        handles.Windows_template{1,1}(end+1,:) = xx;
        handles.Windows_template{2,1}(end+1) = ibeat;
    case 2
        set(h,'facecolor',[.8 0 .8],'EdgeColor',[.8 0 .8]);
        handles.Windows_template{1,2}(end+1,:) = xx;
        handles.Windows_template{2,2}(end+1) = ibeat;
        set(a,'marker','square','markerfacecolor','b')
end

set(handles.axes1,'xlim',[0 size(handles.S,1)]/handles.ParamSig.frequency*1000)
guidata(hObject,handles);


% --- Executes on selection change in tag_select_type.
function tag_select_type_Callback(hObject, eventdata, handles)
% hObject    handle to tag_select_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_select_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_select_type



% --- Executes during object creation, after setting all properties.
function tag_select_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_select_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_ok_button.
function tag_ok_button_Callback(hObject, eventdata, handles)
% hObject    handle to tag_ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% point = str2double(handles.tag_point_number.String);
cc = findobj(0,'name','Inputs for Markers');
data0 = guidata(cc(1));
DeltaT = str2double(data0.tag_DeltaT.String);
% RT min-max
if ~isempty(handles.Windows_template{1,1})
    data0.tag_RTmax.String = num2str(nanmean(handles.Windows_template{1,1}(:,2))-DeltaT,'%1.1f');
    data0.tag_RTmin.String = num2str(nanmean(handles.Windows_template{1,1}(:,1))-DeltaT,'%1.1f');
end

% QRS
if ~isempty(handles.Windows_template{1,2})
    data0.tag_ATmax.String = num2str(nanmean(handles.Windows_template{1,2}(:,2))-DeltaT,'%1.1f');
    data0.tag_DeltaT.String = num2str(DeltaT-nanmean(handles.Windows_template{1,2}(:,1)),'%1.1f');
end



% assignin('base','Template',Template);
close(handles.figure1)

function tag_point_number_Callback(hObject, eventdata, handles)
% hObject    handle to tag_point_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_point_number as text
%        str2double(get(hObject,'String')) returns contents of tag_point_number as a double


% --- Executes during object creation, after setting all properties.
function tag_point_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_point_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tab_plot_button.
function tab_plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to tab_plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ibeat = str2double(handles.tag_point_number.String);
ichan = handles.listbox1.Value;
handles.tag_slider.Value =ibeat;

ax = [handles.axes1];
hold(ax(1),'off');
t = [0:size(handles.S,1)-1]/handles.ParamSig.frequency*1000;
plot(handles.axes1,t,squeeze(handles.S(:,ibeat,ichan)));
xlabel(ax,'Time (ms)');

ff = findobj(handles.axes2,'color','r','linestyle','--');
delete(ff);
hold(handles.axes2,'on');
plot(handles.axes2,handles.spikes(ibeat)*[1 1]/1000,get(handles.axes2,'ylim'),'--r');

cmap = lines;
if handles.tag_do_show_W.Value
    hold(handles.axes1,'on');
    aa = gobjects(size(handles.Windows_template{1,1},1),1);
    for i = 1:size(handles.Windows_template{1,1},1);
        aa(i,1) = plot(handles.axes1,handles.Windows_template{1,1}(i,1)*[1 1],get(handles.axes1,'ylim'),'--');
        aa(i,2) = plot(handles.axes1,handles.Windows_template{1,1}(i,2)*[1 1],get(handles.axes1,'ylim'),'--');
        set(aa(i,:),'color',cmap(i,:));
    end
    
end
set(ax,'xlim',[0 t(end)]);


% --- Executes on slider movement.
function tag_slider_Callback(hObject, eventdata, handles)
% hObject    handle to tag_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


handles.tag_point_number.String = num2str(round(handles.tag_slider.Value));
tab_plot_button_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function tag_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in tag_delete_leads.
function tag_delete_leads_Callback(hObject, eventdata, handles)
% hObject    handle to tag_delete_leads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = {'I','II','III','aVL','aVR','aVF','V1','V2','V3','V4','V5','V6'};
ii_ko = listdlg('ListString',str,'SelectionMode','Multibple','PromptString','Select leads you want to discard from following analysis');

if ~isempty(ii_ko)
    Leads_ko = str(ii_ko);
    handles.Leads_ko = Leads_ko;
end
guidata(hObject,handles);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_do_show_W.
function tag_do_show_W_Callback(hObject, eventdata, handles)
% hObject    handle to tag_do_show_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_do_show_W


% --- Executes on button press in tag_delete_W.
function tag_delete_W_Callback(hObject, eventdata, handles)
% hObject    handle to tag_delete_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode(handles.figure1);
info_struct = getCursorInfo(dcm_obj);
rr = diff(handles.spikes);
if ~isempty(info_struct)
    t0 = info_struct.Position(1);
    rr0 = info_struct.Position(2);
    ibeat = find(rr==rr0 & t0*1000==handles.spikes(1:end-1));
    
    iiko = handles.Windows_template{2,1}==ibeat;
    handles.Windows_template{1,1}(iiko,:) = [];
    handles.Windows_template{2,1}(iiko) = [];
    
    xx = findobj(handles.axes2,'YData',rr0,'XData',t0);
    delete(xx)
end
%
guidata(hObject, handles);
