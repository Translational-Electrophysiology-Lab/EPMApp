function varargout = PlotSignals(varargin)
% PLOTSIGNALS MATLAB code for PlotSignals.fig
%      PLOTSIGNALS, by itself, creates a new PLOTSIGNALS or raises the existing
%      singleton*.
%
%      H = PLOTSIGNALS returns the handle to a new PLOTSIGNALS or the handle to
%      the existing singleton*.
%
%      PLOTSIGNALS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTSIGNALS.M with the given input arguments.
%
%      PLOTSIGNALS('Property','Value',...) creates a new PLOTSIGNALS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlotSignals_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlotSignals_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlotSignals

% Last Modified by GUIDE v2.5 18-Apr-2014 16:24:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PlotSignals_OpeningFcn, ...
    'gui_OutputFcn',  @PlotSignals_OutputFcn, ...
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


% --- Executes just before PlotSignals is made visible.
function PlotSignals_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlotSignals (see VARARGIN)

% Choose default command line output for PlotSignals
handles.output = hObject;

%% ADD DATA
handles.hObject_main = varargin{1};
data = guidata(handles.hObject_main);
handles.ParamSig = data.ParamSig;
if isfield(data,'signals')
    handles.Sig = data.signals;
else
    handles.Sig = data.signals_raw;
end
if isfield(data,'signals_proc')
    handles.Sig_proc = data.signals_proc;
else
    handles.Sig_proc = [];
end
if isfield(data,'spikes')
    handles.spikes = data.spikes;
end
if isfield(data,'MarkersC')
    handles.MarkersC = data.MarkersC;
end

if isfield(data,'LOG')
    handles.LOG=data.LOG;
end
if isfield(handles,'iscut')
    handles.iscut = 0;
end


set(handles.tag_list,'string',handles.ParamSig.Label)
set(handles.tag_ax_sig,'xtick',[],'ytick',[]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PlotSignals wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PlotSignals_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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


% --- Executes on button press in tag_plot.
function tag_plot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ichan = get(handles.tag_list,'value');
str = get(handles.tag_list,'string');
chan_name = str{ichan};
Scales = [1 2 10 inf];
if get(handles.tag_scale,'value')==1
    set(handles.tag_scale,'value',2)
end
iscale = get(handles.tag_scale,'value')-1;
distchan = -[1:length(ichan)]/Scales(iscale);

if get(handles.tag_filt,'value')
    if ~isempty(handles.Sig_proc)
        Sn = handles.Sig_proc(:,ichan)./repmat(max(abs(handles.Sig_proc(:,ichan))),[size(handles.Sig_proc,1),1])+ repmat(distchan,[size(handles.Sig_proc,1),1]);
    end
else
    Sn = handles.Sig(:,ichan)./repmat(max(abs(handles.Sig(:,ichan))),[size(handles.Sig,1),1])+ repmat(distchan,[size(handles.Sig,1),1]);
end


if get(handles.tag_x_units,'value')==1
    x=[1:size(Sn,1)];xlab = 'samples';
elseif get(handles.tag_x_units,'value')==2
    x=[1:size(Sn,1)]/handles.ParamSig.frequency*1000;xlab = 'Time [ms]';
elseif get(handles.tag_x_units,'value')==3
    x=[1:size(Sn,1)]/handles.ParamSig.frequency;;xlab = 'Time [sec]';
elseif get(handles.tag_x_units,'value')==4
    x=[1:size(Sn,1)]/handles.ParamSig.frequency/60;;xlab = 'Time [min]';
end

if ~isfield(handles,'xlab')
    handles.xlab=[];
end

xx = get(handles.tag_ax_sig,'xlim');
hplot=plot(handles.tag_ax_sig,x,Sn);
if ~isequal(xx,[0 1])&isequal(xlab,handles.xlab)
    set(handles.tag_ax_sig,'xlim',xx,'ylim',[min(Sn(:))-eps max(Sn(:))+eps])
else
    axis(handles.tag_ax_sig,'tight')
end
if iscale~=4
    set(handles.tag_ax_sig,'ytick',distchan(end:-1:1),'yticklabel',handles.ParamSig.Label(ichan(end:-1:1)))
end

ll = legend(hplot,handles.ParamSig.Label(ichan));
set(ll,'units','normalized','position',[.92 .2 .075 .04*length(ichan)],'box','off')
xlabel(handles.tag_ax_sig,xlab)
handles.xlab=xlab;
%update
guidata(hObject,handles);

% --- Executes on button press in tag_cut.
function tag_cut_Callback(hObject, eventdata, handles)
% hObject    handle to tag_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = get(handles.tag_x_units,'string');
xx = get(handles.tag_ax_sig,'xlim');
if isequal(str{get(handles.tag_x_units,'value')},'s')
    ini = max(1,round(xx(1)*handles.ParamSig.frequency));
    fin = round(min(xx(2)*handles.ParamSig.frequency,size(handles.Sig,1)));
elseif isequal(str{get(handles.tag_x_units,'value')},'ms')
    ini = max(1,round(xx(1)/1000*handles.ParamSig.frequency));
    fin = round(min(xx(2)/1000*handles.ParamSig.frequency,size(handles.Sig,1)));
elseif isequal(str{get(handles.tag_x_units,'value')},'min')
    ini = max(1,round(xx(1)*60*handles.ParamSig.frequency));
    fin = round(min(xx(2)*60*handles.ParamSig.frequency,size(handles.Sig,1)));
elseif isequal(str{get(handles.tag_x_units,'value')},'samples')
    ini = max(1,round(xx(1)));
    fin = round(min(xx(2),size(handles.Sig,1)));
end

% ini = max(1,round(xx(1)));fin = round(min(xx(2),size(handles.Sig,1)));
Sig = handles.Sig(ini:fin,:);
handles.ParamSig.CropTime = [ini fin];
if isfield(handles,'LOG')
    handles.LOG = [handles.LOG,['- Cropping sig: [',num2str(ini),':',num2str(fin),']']];
end
%

if isfield(handles,'Sig_proc')
    if ~isempty(handles.Sig_proc)
        handles.Sig_proc = handles.Sig_proc(ini:fin,:);
    end
end

if isfield(handles,'MarkersC')
    vv = fieldnames(handles.MarkersC);
    hh =handles.spikes<ini/handles.ParamSig.frequency*1000 | handles.spikes>fin/handles.ParamSig.frequency*1000;
    for j = 1:length(vv)
        if size(handles.MarkersC.(vv{j}),1)==length(handles.spikes)
            handles.MarkersC.(vv{j})(hh,:) =[];
            handles.MarkersC.(vv{j}) = handles.MarkersC.(vv{j}) - (ini/handles.ParamSig.frequency*1000);
        end
    end
end

if isfield(handles,'spikes')
    handles.spikes(handles.spikes<ini/handles.ParamSig.frequency*1000 | handles.spikes>fin/handles.ParamSig.frequency*1000)=[];
    handles.spikes =handles.spikes - (ini/handles.ParamSig.frequency*1000);
end
set(handles.tag_cut,'backgroundcolor',[.6 .6 .6],'foregroundcolor',[1 0 0])
pause(0.25)
set(handles.tag_cut,'backgroundcolor',[1 1 1]*.941,'foregroundcolor',[0 0 0])
% update
handles.iscut =1;
handles.Sig=Sig;
handles.xlab = []; % just to resize
guidata(hObject, handles);

% --- Executes on button press in tag_done.
function tag_done_Callback(hObject, eventdata, handles)
% hObject    handle to tag_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = guidata(handles.hObject_main);
if isfield(data,'LOG')
    data.LOG = [data.LOG,['- Visualization ended']];
end

if handles.iscut
    data.signals = handles.Sig;
    data.signals_proc = handles.Sig_proc;
    if isfield(handles,'MarkersC')
        data.MarkersC = handles.MarkersC;
    end
    
    data.ParamSig = handles.ParamSig;
    if isfield(data.ParamSig,'Times') % eliminate useless data from ensite
        data.ParamSig=rmfield(data.ParamSig,'Times');
    end
    if isfield(data.ParamSig,'Times_name') % eliminate useless data from ensite
        data.ParamSig=rmfield(data.ParamSig,'Times_name');
    end
    if isfield(handles.ParamSig,'CropTime')&isfield(data,'LOG')
        data.LOG = [data.LOG,['- Signals cropped : time [',num2str(handles.ParamSig.CropTime),']']];
    end
    if isfield(handles.ParamSig,'CropChan')&isfield(data,'LOG')
        data.LOG = [data.LOG,['- Signals cropped : saved channels [',num2str(handles.ParamSig.CropChan),']']];
    end
    if get(data.tag_save_session,'value')&~isempty(data.outputname)
        signals = handles.Sig;
        save(data.outputname,'signals','-append')
        %
        if isfield(data,'LOG')
            data.LOG = [data.LOG,['- saved : signals_raw']];
        end
    end
    if isfield(handles,'spikes')
        data.spikes = handles.spikes;
    end
end
if isfield(data,'LOG')
    set(data.tag_Log,'string',data.LOG);
end
% update
guidata(handles.hObject_main,data);
close(handles.figure1)

% --- Executes on selection change in tag_scale.
function tag_scale_Callback(hObject, eventdata, handles)
% hObject    handle to tag_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_scale contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_scale


% --- Executes during object creation, after setting all properties.
function tag_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tag_x_units.
function tag_x_units_Callback(hObject, eventdata, handles)
% hObject    handle to tag_x_units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_x_units contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_x_units


% --- Executes during object creation, after setting all properties.
function tag_x_units_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_x_units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_filt.
function tag_filt_Callback(hObject, eventdata, handles)
% hObject    handle to tag_filt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_filt


% --- Executes on button press in tag_export.
function tag_export_Callback(hObject, eventdata, handles)
% hObject    handle to tag_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fig2 = figure;
copyobj(handles.tag_ax_sig, Fig2);


% --- Executes on button press in tag_cut_channels.
function tag_cut_channels_Callback(hObject, eventdata, handles)
% hObject    handle to tag_cut_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iiko = get(handles.tag_list,'value');
str = get(handles.tag_list,'string');
iiok = [1:length(str)];iiok(iiko)=[];

str2 = cell(length(str),1);
for i = 1:length(str)
    if ~isempty(find(iiko==i))
        str2(i) = {[str{i},' (cancelled)']};
    else
        str2(i) = {[str{i},' (ok)']};
    end
end
[answ,ok] = listdlg('liststring',str2,'promptstring','Do you want to cancel selected channels?','listsize',[310 500],'initialvalue',iiko);
if ok
    handles.Sig(:,iiko)=[];
    if ~isempty(handles.Sig_proc)
        handles.Sig_proc(:,iiko)=[];
    end
    handles.ParamSig.Label(iiko)=[];
    set(handles.tag_list,'string',handles.ParamSig.Label,'value',[1:length(handles.ParamSig.Label)])
    handles.ParamSig.CropChan = [iiok];
    handles.iscut=true;
end
% update
guidata(hObject,handles)

