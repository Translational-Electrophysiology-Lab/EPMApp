function varargout = GUI_plot_carto_point(varargin)
% GUI_PLOT_CARTO_POINT MATLAB code for GUI_plot_carto_point.fig
%      GUI_PLOT_CARTO_POINT, by itself, creates a new GUI_PLOT_CARTO_POINT or raises the existing
%      singleton*.
%
%      H = GUI_PLOT_CARTO_POINT returns the handle to a new GUI_PLOT_CARTO_POINT or the handle to
%      the existing singleton*.
%
%      GUI_PLOT_CARTO_POINT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PLOT_CARTO_POINT.M with the given input arguments.
%
%      GUI_PLOT_CARTO_POINT('Property','Value',...) creates a new GUI_PLOT_CARTO_POINT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_plot_carto_point_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_plot_carto_point_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_plot_carto_point

% Last Modified by GUIDE v2.5 19-Jun-2017 21:29:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_plot_carto_point_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_plot_carto_point_OutputFcn, ...
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


% --- Executes just before GUI_plot_carto_point is made visible.
function GUI_plot_carto_point_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_plot_carto_point (see VARARGIN)

if nargin>0
    dir_load = varargin{1};
else
    dir_load = uigetdir(pwd,'Select a folder containing Carto Data');
    
end
dir_load = [dir_load,'\'];
% Choose default command line output for GUI_plot_carto_point
handles.output = hObject;
% ===
X = dir([dir_load,'*ECG_Export*.txt']);
Map = cell(1,length(X));
for i = 1:length(X)
    ii = find(X(i).name=='_');
    Map{i} = X(i).name(1:ii(1)-1);
end
Map_type = unique(Map);
% -
istot = [];
is_name_maps = cell(1,length(Map_type));
points_name_maps = cell(1,length(Map_type));
for jmap = 1:length(Map_type)
    imap  = find(~cellfun(@isempty,strfind({X.name},Map_type{jmap})));
    n = nan(1,length(imap));
    for i = 1:length(imap)
        ii = find(X(imap(i)).name=='_');
        n(i) = str2double(X(imap(i)).name(ii(1)+2:ii(2)-1));
    end
    [~,is] = sort(n);
    istot = [istot imap(is)];
    is_name_maps(jmap) = {imap(is)};
    points_name_maps{jmap} = {X(imap(is)).name};
end
points_name = {X(istot).name};
clear n imap jmap
%%
handles.tag_list_points.String = points_name;
handles.tag_list_points.Max = length(points_name);
handles.mat.dir_load = dir_load;
handles.mat.points_name = points_name;
handles.mat.points_name_maps = points_name;
handles.mat.Map_type = Map_type;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_plot_carto_point wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_plot_carto_point_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in tag_list_points.
function tag_list_points_Callback(hObject, eventdata, handles)
% hObject    handle to tag_list_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_list_points contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_list_points


% --- Executes during object creation, after setting all properties.
function tag_list_points_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_list_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_plot.
function tag_plot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


xxl = get(handles.tag_ax,'xlim');
if diff(xxl)==1
    xxl = [0 2.5];
end
yyl = get(handles.tag_ax,'ylim');
do_original = handles.tag_Yaxis.Value==4;
ii = handles.tag_list_points.Value;
% points
leg_point = handles.tag_list_points.String(ii);
for j = 1:length(leg_point)
    jj =find(leg_point{j}=='_');
    leg_point{j} =leg_point{j}((jj(1)+1):(jj(2)-1));
end
clear j jj

sig = cell(1,length(ii));
for jp = 1:length(ii)
    loadname = [handles.mat.dir_load,handles.tag_list_points.String{ii(jp)}];
    [s,labels] = LoadWaveformsFromCarto(loadname);
    sig{jp} = s;
end
clear s

[cc,~,iichan] = intersect(handles.tag_list_chan.String(handles.tag_list_chan.Value),labels,'stable');

splot = nan(size(sig{1},1),length(iichan)*length(ii));
leg = cell(1,length(iichan)*length(ii));
ll = cell(1,length(iichan));
for jp = 1:length(ii)
    splot(:,[1:length(iichan)]+(jp-1)*length(iichan)) = sig{jp}(:,iichan);
    for jp2 = 1:length(iichan)
        a = cc{jp2}(1 : find(cc{jp2}=='(')-1);
        ll(jp2) = {[leg_point{jp},'-',a]};
        ll{jp2}(ll{jp2}=='_') = '-';
    end
    leg([1:length(iichan)]+(jp-1)*length(iichan)) = ll;
end
% reorder splot
jj = [reshape([1:length(ii)*length(iichan)],[length(iichan) length(ii)])]';
jj = jj(:);

leg = leg(jj);
splot = splot(:,jj);
clear jj

if ~do_original
    splot = splot./(ones(size(splot,1),1)*mad(splot,1))/3 + ones(size(splot,1),1)*[0:-1:-(size(splot,2)-1)];
    handles.tag_Yaxis.Value = 1;
end

plot(handles.tag_ax,[1:size(splot,1)]/1000,splot);
legend(leg)
set(handles.tag_ax,'box','off')
set(handles.tag_ax,'xlim',xxl)
% set(handles.tag_ax,'ylim',yyl)





% --- Executes on selection change in tag_list_chan.
function tag_list_chan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_list_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_list_chan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_list_chan


% --- Executes during object creation, after setting all properties.
function tag_list_chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_list_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_select_P.
function tag_select_P_Callback(hObject, eventdata, handles)
% hObject    handle to tag_select_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



ii = handles.tag_list_points.Value;
Ltot = cell(1,length(ii));
for jp = 1:length(ii)
    loadname = [handles.mat.dir_load,handles.tag_list_points.String{ii(jp)}];
    [~,L,I] = LoadWaveformsFromCarto(loadname);
    Ltot{jp} = L(:);
end
control_points = zeros(1,length(jp));
for jp = 1:length(ii)
    control_points(jp) = isequal(Ltot{jp},Ltot{1});
end
if sum(control_points~=1)>0
    error('Points must contain recordings from the same channels');
end
clear jp
load ListOfChannels;

[cc,ia] = intersect(ListOfChannels_Template,Ltot{1},'stable');

handles.tag_list_chan.String = [];
pause(0.2)
handles.tag_list_chan.String = cc;
handles.tag_list_chan.Max = length(cc);
handles.mat.chan_name = cc;
handles.mat.chan_order = ia;
% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in tag_Yaxis.
function tag_Yaxis_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Yaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_Yaxis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_Yaxis
zoom('off')
dataObjs = get(handles.tag_ax, 'Children');
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = [get(dataObjs, 'YData')];
ydata = ydata(end:-1:1);
y = cell2mat(ydata).';
x = cell2mat(xdata).';
% g = round(nanmean(y));
d = (mean((diff(nanmean(y)))));
g = [-size(y,2)*d:d:-d]+1;
ll = legend;
ls = ll.String;
if handles.tag_Yaxis.Value==2
    xx = get(handles.tag_ax,'xlim');
    y2 = y + ones(size(y,1),1)*(g(:).')/2;
    hold(handles.tag_ax,'off');
    plot(handles.tag_ax,x(:,1),y2);
    set(handles.tag_ax,'xlim',xx);
    set(handles.tag_ax,'ylim',[min(y2(:)) max(y2(:))]);
elseif handles.tag_Yaxis.Value==3
    xx = get(handles.tag_ax,'xlim');
    y2 = y - ones(size(y,1),1)*(g(:).')/2;
    hold(handles.tag_ax,'off');
    plot(handles.tag_ax,x(:,1),y2);
    set(handles.tag_ax,'xlim',xx);
    set(handles.tag_ax,'ylim',[min(y2(:)) max(y2(:))]);
    
elseif  handles.tag_Yaxis.Value==1|handles.tag_Yaxis.Value==4
    tag_plot_Callback(hObject,[],handles);
elseif  handles.tag_Yaxis.Value==5
    xx = get(handles.tag_ax,'xlim');
    yy = get(handles.tag_ax,'ylim');
    axis(handles.tag_ax,'tight');
    zoom(handles.tag_ax,'reset')
    set(handles.tag_ax,'ylim',yy);
    set(handles.tag_ax,'xlim',xx);
end

% axis tight
set(handles.tag_ax,'box','off')
legend(handles.tag_ax,ls)

% --- Executes during object creation, after setting all properties.
function tag_Yaxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_Yaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
