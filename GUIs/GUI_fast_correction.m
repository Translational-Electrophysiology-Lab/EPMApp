function varargout = GUI_fast_correction(varargin)
% GUI_FAST_CORRECTION MATLAB code for GUI_fast_correction.fig
%      GUI_FAST_CORRECTION, by itself, creates a new GUI_FAST_CORRECTION or raises the existing
%      singleton*.
%
%      H = GUI_FAST_CORRECTION returns the handle to a new GUI_FAST_CORRECTION or the handle to
%      the existing singleton*.
%
%      GUI_FAST_CORRECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FAST_CORRECTION.M with the given input arguments.
%
%      GUI_FAST_CORRECTION('Property','Value',...) creates a new GUI_FAST_CORRECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_fast_correction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_fast_correction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_fast_correction

% Last Modified by GUIDE v2.5 19-Jun-2017 17:48:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_fast_correction_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_fast_correction_OutputFcn, ...
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


% --- Executes just before GUI_fast_correction is made visible.
function GUI_fast_correction_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_fast_correction (see VARARGIN)

% Choose default command line output for GUI_fast_correction
handles.output = hObject;

%% ADD DATA
data = varargin{1};
handles.mat.spikes = data.spikes(:);
handles.mat.S = data.signals;
if isfield(data,'signals_proc')
    handles.mat.signals_proc = data.signals_proc;
end
if isfield(data,'SNR')
    handles.mat.SNR = data.SNR;
end
if isfield(data,'Markers')
    handles.mat.Markers = data.Markers;
end
if isfield(data,'MarkersC')
    handles.mat.MarkersC = data.MarkersC;
end
if isfield(data,'sig_corr')
    handles.mat.sig_corr = data.sig_corr;
end
if isfield(data,'sig_corr_beats')
    handles.mat.sig_corr_beats = data.sig_corr_beats;
end

handles.mat.ParamSig = data.ParamSig;

List = handles.mat.ParamSig.Label;
for k = 1:length(List)
    List{k} = ['[',num2str(k),']',List{k}];
end
if isfield(handles.mat,'SNR')
    for j=1:length(handles.mat.SNR)
        List{j} = [List{j},'(',num2str(handles.mat.SNR(j),2),'dB)'];
    end
end
if isfield(handles.mat,'sig_corr')
    if ~isempty(handles.mat.sig_corr)
        for j=1:length(List)
            List{j} = [List{j},'/[',num2str(handles.mat.sig_corr(j),2),']'];
        end
    end
end
set(handles.tag_list,'string',List)
set(handles.tag_ax,'xticklabel',[],'yticklabel',[])
set(handles.tag_Nchan,'string','10')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_fast_correction wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_fast_correction_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in tag_marker_type.
function tag_marker_type_Callback(hObject, eventdata, handles)
% hObject    handle to tag_marker_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_marker_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_marker_type


% --- Executes during object creation, after setting all properties.
function tag_marker_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_marker_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_plot.
function tag_plot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hold(handles.tag_ax,'off');
iic = handles.tag_list.Value;
splot = handles.mat.signals_proc(:,iic);
splot = splot./(ones(size(splot,1),1)*mad(splot,1))/3 + ones(size(splot,1),1)*[0:-1:-(size(splot,2)-1)];
handles.tag_Yaxis.Value = 1;

t = [1:size(splot,1)]/handles.mat.ParamSig.frequency*1000;
aa = plot(handles.tag_ax,t,splot);
xlabel(handles.tag_ax,'Time (ms)');

leg = cell(1,length(iic));
for j=1:length(iic)
    leg{j} = num2str(iic(j));
end
legend(aa,leg)
set(handles.tag_ax,'ytick',[-size(splot,2)+1:1:0],'yticklabel',leg(end:-1:1))
% - plot markers
if handles.tag_marker_type.Value==1
    x = handles.mat.MarkersC.dt(:,iic);
elseif handles.tag_marker_type.Value==2
    x = handles.mat.MarkersC.rt_Wyatt(:,iic);
end
hold(handles.tag_ax,'on');
ix = round(x/handles.mat.ParamSig.frequency*1000);
for j = 1:length(iic)
    ip = ix(~isnan(ix(:,j)),j);
    plot(handles.tag_ax,t(ip),splot(ip,j),'x','color',aa(j).Color,'markersize',8,'linewidth',2)
end


% --- Executes on slider movement.
function tag_slider_Callback(hObject, eventdata, handles)
% hObject    handle to tag_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function tag_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function tag_Nchan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Nchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_Nchan as text
%        str2double(get(hObject,'String')) returns contents of tag_Nchan as a double


% --- Executes during object creation, after setting all properties.
function tag_Nchan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_Nchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tag_list.
function tag_list_Callback(hObject, eventdata, handles)
% hObject    handle to tag_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_list


% --- Executes during object creation, after setting all properties.
function tag_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_apply_Nchan.
function tag_apply_Nchan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_apply_Nchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in tag_Yaxis.
function tag_Yaxis_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Yaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_Yaxis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_Yaxis
zoom('off')
dataObjs = get(handles.tag_ax, 'Children');
iiobj = findobj(dataObjs,'Marker','none');
xdata = get(iiobj, 'XData');  %data from low-level grahics objects
ydata = [get(iiobj, 'YData')];
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
    aa = plot(handles.tag_ax,x(:,1),y2);
    set(handles.tag_ax,'xlim',xx);
    set(handles.tag_ax,'ylim',[min(y2(:)) max(y2(:))]);
elseif handles.tag_Yaxis.Value==3
    xx = get(handles.tag_ax,'xlim');
    y2 = y - ones(size(y,1),1)*(g(:).')/2;
    hold(handles.tag_ax,'off');
    aa =plot(handles.tag_ax,x(:,1),y2);
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

% - plot markers
t = [1:size(y2,1)]/handles.mat.ParamSig.frequency*1000;
iic = handles.tag_list.Value;
if handles.tag_marker_type.Value==1
    x = handles.mat.MarkersC.dt(:,iic);
elseif handles.tag_marker_type.Value==2
    x = handles.mat.MarkersC.rt_Wyatt(:,iic);
end
hold(handles.tag_ax,'on');
ix = round(x/handles.mat.ParamSig.frequency*1000);
for j = 1:length(iic)
    ip = ix(~isnan(ix(:,j)),j);
    plot(handles.tag_ax,t(ip),y2(ip,j),'x','color',aa(j).Color,'markersize',8,'linewidth',2)
end
set(handles.tag_ax,'ytick',[-size(y2,2)+1:1:0],'yticklabel',ls(end:-1:1))

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



function tag_win_width_Callback(hObject, eventdata, handles)
% hObject    handle to tag_win_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_win_width as text
%        str2double(get(hObject,'String')) returns contents of tag_win_width as a double


% --- Executes during object creation, after setting all properties.
function tag_win_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_win_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_modify.
function tag_modify_Callback(hObject, eventdata, handles)
% hObject    handle to tag_modify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in tag_all_beats.
function tag_all_beats_Callback(hObject, eventdata, handles)
% hObject    handle to tag_all_beats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_all_beats


% --- Executes on button press in tag_all_chan.
function tag_all_chan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_all_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_all_chan


% --- Executes on button press in tag_continuous_click.
function tag_continuous_click_Callback(hObject, eventdata, handles)
% hObject    handle to tag_continuous_click (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_continuous_click
