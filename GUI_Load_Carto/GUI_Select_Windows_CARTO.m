function varargout = GUI_Select_Windows_CARTO(varargin)
% GUI_SELECT_WINDOWS_CARTO MATLAB code for GUI_Select_Windows_CARTO.fig
%      GUI_SELECT_WINDOWS_CARTO, by itself, creates a new GUI_SELECT_WINDOWS_CARTO or raises the existing
%      singleton*.
%
%      H = GUI_SELECT_WINDOWS_CARTO returns the handle to a new GUI_SELECT_WINDOWS_CARTO or the handle to
%      the existing singleton*.
%
%      GUI_SELECT_WINDOWS_CARTO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SELECT_WINDOWS_CARTO.M with the given input arguments.
%
%      GUI_SELECT_WINDOWS_CARTO('Property','Value',...) creates a new GUI_SELECT_WINDOWS_CARTO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Select_Windows_CARTO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Select_Windows_CARTO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Select_Windows_CARTO

% Last Modified by GUIDE v2.5 18-May-2019 14:34:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Select_Windows_CARTO_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Select_Windows_CARTO_OutputFcn, ...
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


% --- Executes just before GUI_Select_Windows_CARTO is made visible.
function GUI_Select_Windows_CARTO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Select_Windows_CARTO (see VARARGIN)

% Choose default command line output for GUI_Select_Windows_CARTO
handles.output = hObject;
handles.S = varargin{1};
handles.Windows_template = cell(1,4);

handles.tag_slider.SliderStep = [1 10]/size(handles.S,3);
handles.tag_slider.Max = size(handles.S,3);
handles.tag_slider.Min = 1;
handles.tag_slider.Value = 1;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes GUI_Select_Windows_CARTO wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Select_Windows_CARTO_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = [];

% --- Executes on button press in tag_select_button.
function tag_select_button_Callback(hObject, eventdata, handles)
% hObject    handle to tag_select_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ax = [handles.axes1 handles.axes2];

xx = get(ax(1),'xlim');
xx(1) = max([0,xx(1)]);
xx(2) = min([xx(2) size(handles.S,1)]);

hold(ax(1),'on')
h(1) = patch([xx(1) xx(2) xx(2) xx(1)],[-6 -6 1 1],[0 .8 .8],'parent',ax(1));
h(2) = patch([xx(1) xx(2) xx(2) xx(1)],[-12 -12 -5 -5],[0 .8 .8],'parent',ax(2));
set(h,'FaceAlpha',0.1);
set(h,'EdgeColor','k');
set(h,'LineWidth',2);

% 
is = handles.tag_select_type.Value;
if is==1
    set(h,'facecolor',[0 .8 .8],'EdgeColor',[0 .8 .8]);
     handles.Windows_template{1} = xx;
elseif is==2
    set(h,'facecolor',[.8 0 .8],'EdgeColor',[.8 0 .8]);
    handles.Windows_template{2} = xx;
elseif is==3
    set(h,'facecolor',[.8 .1 .1],'EdgeColor',[.8 0 .8]);
    handles.Windows_template{3} = xx;
elseif is==4
    set(h,'facecolor',[0 .8 0],'EdgeColor',[0 .8 0]);
    handles.Windows_template{4} = xx;
end
    
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

point = str2double(handles.tag_point_number.String);
data = cell(size(handles.Windows_template));
for i = 1:4
    if ~isempty(handles.Windows_template{i})
data{i} = [round(handles.Windows_template{i}(1)) : round(handles.Windows_template{i}(2))];  
    end
end
Template.point = point;
Template.Windows = data;
if isfield(handles,'Leads_ko')
    Template.Lead_ko = handles.Leads_ko;
end

% assignin('base','Template',Template);

cc = findobj(0,'name','GUI_Load_Carto')
hmain = guidata(cc);
hmain.mat.Template = Template;
guidata(hmain.figure1,hmain);

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

ipoint = str2double(handles.tag_point_number.String);
handles.tag_slider.Value =ipoint;

ax = [handles.axes1 handles.axes2];
hold(ax(1),'off')
hold(ax(2),'off')
plot(handles.axes1,handles.S(:,1:6,ipoint));
plot(handles.axes2,handles.S(:,7:12,ipoint));
set(ax,'xlim',[0 2500]);
set(ax(1),'ytick',[-5:0],'yticklabel',{'aVF','aVR','aVL','III','II','I'});
set(ax(2),'ytick',[-11:-6],'yticklabel',{'V6','V5','V4','V3','V2','V1'});
linkaxes(ax,'x')

% 
colors = [0 .8 .8;.8 0 .8;.8 .1 .1;0 .8 0];
for is = 1:4;
    if ~isempty(handles.Windows_template{is})
        xx = handles.Windows_template{is};
        hold(ax(1),'on');hold(ax(2),'on')
        h(1) = patch([xx(1) xx(2) xx(2) xx(1)],[-6 -6 1 1],[0 .8 .8],'parent',ax(1));
        h(2) = patch([xx(1) xx(2) xx(2) xx(1)],[-12 -12 -5 -5],[0 .8 .8],'parent',ax(2));
        set(h,'FaceAlpha',0.1);
        set(h,'EdgeColor',colors(is,:),'facecolor',colors(is,:));
        set(h,'LineWidth',2);
        
    end
end

    


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


% --------------------------------------------------------------------
function uitoggletool1_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom off
zoom on


% --------------------------------------------------------------------
function uitoggletool2_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom off
zoom out
