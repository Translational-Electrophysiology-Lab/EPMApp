function varargout = Check_Channels_GUI_ATRT(varargin)
% CHECK_CHANNELS_GUI_ATRT MATLAB code for Check_Channels_GUI_ATRT.fig
%      CHECK_CHANNELS_GUI_ATRT, by itself, creates a new CHECK_CHANNELS_GUI_ATRT or raises the existing
%      singleton*.
%
%      H = CHECK_CHANNELS_GUI_ATRT returns the handle to a new CHECK_CHANNELS_GUI_ATRT or the handle to
%      the existing singleton*.
%
%      CHECK_CHANNELS_GUI_ATRT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHECK_CHANNELS_GUI_ATRT.M with the given input arguments.
%
%      CHECK_CHANNELS_GUI_ATRT('Property','Value',...) creates a new CHECK_CHANNELS_GUI_ATRT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Check_Channels_GUI_ATRT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Check_Channels_GUI_ATRT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Check_Channels_GUI_ATRT

% Last Modified by GUIDE v2.5 01-Nov-2017 12:22:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Check_Channels_GUI_ATRT_OpeningFcn, ...
    'gui_OutputFcn',  @Check_Channels_GUI_ATRT_OutputFcn, ...
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


% --- Executes just before Check_Channels_GUI_ATRT is made visible.
function Check_Channels_GUI_ATRT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Check_Channels_GUI_ATRT (see VARARGIN)

% Choose default command line output for Check_Channels_GUI_ATRT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Check_Channels_GUI_ATRT wait for user response (see UIRESUME)
% uiwait(handles.figure1);
data = varargin{1};



%% = Parameters
if isfield(data.ParamSig,'chan_deleted')
    data.chan_deleted = data.ParamSig.chan_deleted;
    data.ParamPreProc = data.ParamSig.chan_deleted.ParamPreProc;
else
    data.ParamPreProc.SNR_th = 0;
    data.ParamPreProc.sig_corr_th = 0.75;
    data.ParamPreProc.amp_th = 20;
    data.ParamPreProc.cv_thr = 0.3;

    data.chan_deleted.Amplitude = mean(abs(detrend(data.S,'constant'))>data.ParamPreProc.amp_th).'>0.05;
    data.chan_deleted.SNR = data.SNR(:) < data.ParamPreProc.SNR_th;
    data.chan_deleted.Sig_Corr = mean(data.sig_corr<data.ParamPreProc.sig_corr_th,2)>0.66;
    data.chan_deleted.Deleted_manually = [data.chan_deleted.Amplitude | data.chan_deleted.SNR | data.chan_deleted.Sig_Corr];
    data.chan_deleted.Checked_manually  = zeros(size(data.S,2),1);
    if isfield(data,'Markers')&&isfield(data,'spikes')
        dt = data.Markers.dt - data.spikes(:)*ones(1,size(data.Markers.dt,2));
        rt = data.Markers.rt_Wyatt - data.spikes(:)*ones(1,size(data.Markers.rt_Wyatt,2));
        data.chan_deleted.Checked_manually = mad(rt,1)*1.48./nanmedian(rt,1)>data.ParamPreProc.cv_thr | mad(dt,1)*1.48./nanmedian(dt,1)>data.ParamPreProc.cv_thr;
    end
end
% Table
dat =  cell(size(data.SNR(:),1),6);
dat(:,2) = num2cell([1:size(dat,1)]');
dat(:,1) = data.ParamSig.Label(1:size(dat,1));
if isfield(data,'SNR')
    dat(:,3) = num2cell(data.SNR(:));
end
if isfield(data,'sig_corr')
    dat(:,4) = num2cell(nanmean(data.sig_corr,2));
end
dat(:,5) = num2cell(logical(data.chan_deleted.Deleted_manually));
dat(:,6) = num2cell(logical(data.chan_deleted.Checked_manually));

%
dat = [{'#','Chan','SNR','Corr','Delete','Check'};dat];

set(handles.tag_table,'data',dat,'columneditable',[false false false false true true]);
handles.dat = dat;
set(handles.tag_n_bad_chan,'string',num2str(sum(data.chan_deleted.Deleted_manually)));
set(handles.tag_n_ok_chan,'string',num2str(sum(~data.chan_deleted.Deleted_manually)));
set(handles.tag_n_check,'string',num2str(sum(data.chan_deleted.Checked_manually)));

if isfield(data,'ParamPreProc')
    if isfield(data.ParamPreProc,'SNR_th')
        set(handles.tag_SNR_th,'string',num2str(data.ParamPreProc.SNR_th,2));
    end
    if isfield(data.ParamPreProc,'sig_corr_th')
        set(handles.tag_sig_corr_th,'string',num2str(data.ParamPreProc.sig_corr_th,2));
    end
    if isfield(data.ParamPreProc,'amp_th')
        set(handles.tag_Amp_th,'string',num2str(data.ParamPreProc.amp_th,2));
    end
    if isfield(data.ParamPreProc,'cv_thr')
        set(handles.tag_cv_thr,'string',num2str(data.ParamPreProc.cv_thr,2));
    end
end



handles.data = data;
guidata(hObject,handles);

set(handles.tag_table,'CellSelectionCallback',{@myfun_CellSelectionCallback_ATRT})

% --- Outputs from this function are returned to the command line.
function varargout = Check_Channels_GUI_ATRT_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on selection change in tag_sort_type.
function tag_sort_type_Callback(hObject, eventdata, handles)
% hObject    handle to tag_sort_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_sort_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_sort_type


% --- Executes during object creation, after setting all properties.
function tag_sort_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_sort_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_sort.
function tag_sort_Callback(hObject, eventdata, handles)
% hObject    handle to tag_sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dat_ori = handles.tag_table.Data;
dat = cell(size(dat_ori));
dat(1,:) = dat_ori(1,:);
sort_type = get(handles.tag_sort_type,'string');
% Name
if isequal(sort_type(get(handles.tag_sort_type,'value'),1:4),'Name')
    [~,is] = sort(cell2mat(dat_ori(2:end,2)));
    dat(2:end,:) = dat_ori(is+1,:);
end

% Delete
if isequal(sort_type(get(handles.tag_sort_type,'value'),1:6),'Delete')
    ii1 = find(cell2mat(dat_ori(2:end,5)));
    ii2 = find(~cell2mat(dat_ori(2:end,5)));
    dat(2:length(ii1)+1,:) = dat_ori(ii1+1,:);
    dat(length(ii1)+2:end,:) = dat_ori(ii2+1,:);
end

% Missing
if isequal(sort_type(get(handles.tag_sort_type,'value'),1:7),'Missing')
    ii1 = find(cell2mat(dat_ori(2:end,6)));
    ii2 = find(~cell2mat(dat_ori(2:end,6)));
    dat(2:length(ii1)+1,:) = dat_ori(ii1+1,:);
    dat(length(ii1)+2:end,:) = dat_ori(ii2+1,:);
end

% SNR
if isequal(sort_type(get(handles.tag_sort_type,'value'),1:3),'SNR')
    [~,is] = sort(cell2mat(dat_ori(2:end,3)),'descend');
    dat(2:end,:) = dat_ori(is+1,:);
end
% Corr
if isequal(sort_type(get(handles.tag_sort_type,'value'),1:4),'Corr')
    [~,is] = sort(cell2mat(dat_ori(2:end,4)),'descend');
    dat(2:end,:) = dat_ori(is+1,:);
end
set(handles.tag_table,'data',dat,'columneditable',[false false false true true]);

handles.dat = dat;
guidata(hObject, handles);


% --- Executes on button press in tag_update.
function tag_update_Callback(hObject, eventdata, handles)
% hObject    handle to tag_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.tag_update,'backgroundColor',[1 0.9 0.4]);
pause(0.5);

handles.dat = get(handles.tag_table,'data');
% [dum,iis] = sort(cell2mat(dat(2:end,2)));
% dat = [dat(1,:);dat(iis+1,:)];
iind =  cell2mat(handles.dat(2:end,2));
ii_chan_ko = cell2mat(handles.dat(2:end,5));
iiko = iind(logical(ii_chan_ko));

% Update total number of bad channels
iiko_logic = zeros(size(handles.dat,1)-1,1);
iiko_logic(iiko) = 1;
handles.data.chan_deleted.Deleted_manually = iiko_logic;

% Update total number of channels to be checked
ii_chan_check = cell2mat(handles.dat(2:end,6));
iicheck = iind(logical(ii_chan_check));
iicheck_logic = zeros(size(handles.dat,1)-1,1);
iicheck_logic(iicheck) = 1;
handles.data.chan_deleted.Checked_manually = iicheck_logic;


% transfer it to main
dataori = guidata(handles.data.handle_ori);
handles.data.chan_deleted.ParamPreProc = handles.data.ParamPreProc;
dataori.ParamSig.chan_deleted = handles.data.chan_deleted;

% delete markers
vv = fieldnames(dataori.MarkersC);
iiko = logical(dataori.ParamSig.chan_deleted.Deleted_manually);
for iv = 1:length(vv)
    if size(dataori.MarkersC.(vv{iv}),2)==length(iiko)&size(dataori.MarkersC.(vv{iv}),1)==length(dataori.spikes)
        dataori.MarkersC.(vv{iv})(:,iiko)=nan;
    end
end

%
%%
cc = findobj(0,'name','Markers_GUI_v3');
handles_markers = guidata(cc(1));
L = handles_markers.tag_List.String;
for i = 1:length(L)
    ii1 = find(L{i}=='<');
    ii2 = find(L{i}=='>');
    if ~isempty(ii1)&&~isempty(ii2)
        L{(i)}([1:ii2(2) ii1(3):end]) = [];
    end
end
% =
clear ii*
ii1 = handles.data.chan_deleted.Checked_manually;
ii2 = handles.data.chan_deleted.Deleted_manually;

ii_del = find(ii2);
ii_check = find(logical(ii1)&~logical(ii2));

for i = 1:length(ii_check)
    L{ii_check(i)} = ['<HTML><BODY bgcolor="orange">',L{ii_check(i)},'</BODY></HTML>'];
end

for i = 1:length(ii_del)
    L{ii_del(i)} = ['<HTML><BODY bgcolor="red">',L{ii_del(i)},'</BODY></HTML>'];
end

handles_markers.tag_List.String = L;


% update GUI_SBM_v2.m (field "Manually_checked" added/modified in chan_deleted)
guidata(handles.data.handle_ori,dataori);

cc = findobj(0,'name','check channels');
if ~isempty(cc)
    close(cc)
end
%%
set(handles.tag_update,'backgroundColor',[1 1 1]*.94);




% --- Executes on button press in tag_click_plot.
function tag_click_plot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_click_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_click_plot

if get(handles.tag_click_plot,'value')==1
    % This enables to select a channel and plot in the other window
    set(handles.tag_table,'CellSelectionCallback',{@myfun_CellSelectionCallback_ATRT})
    
else
    % This enables to select a channel and plot in the other window
    set(handles.tag_table,'CellSelectionCallback',[])
    
    
end



function tag_n_bad_chan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_n_bad_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_n_bad_chan as text
%        str2double(get(hObject,'String')) returns contents of tag_n_bad_chan as a double


% --- Executes during object creation, after setting all properties.
function tag_n_bad_chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_n_bad_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_n_ok_chan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_n_ok_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_n_ok_chan as text
%        str2double(get(hObject,'String')) returns contents of tag_n_ok_chan as a double


% --- Executes during object creation, after setting all properties.
function tag_n_ok_chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_n_ok_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_n_check_Callback(hObject, eventdata, handles)
% hObject    handle to tag_n_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_n_check as text
%        str2double(get(hObject,'String')) returns contents of tag_n_check as a double


% --- Executes during object creation, after setting all properties.
function tag_n_check_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_n_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_SNR_th_Callback(hObject, eventdata, handles)
% hObject    handle to tag_SNR_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_SNR_th as text
%        str2double(get(hObject,'String')) returns contents of tag_SNR_th as a double


% --- Executes during object creation, after setting all properties.
function tag_SNR_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_SNR_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_sig_corr_th_Callback(hObject, eventdata, handles)
% hObject    handle to tag_sig_corr_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_sig_corr_th as text
%        str2double(get(hObject,'String')) returns contents of tag_sig_corr_th as a double


% --- Executes during object creation, after setting all properties.
function tag_sig_corr_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_sig_corr_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_Amp_th_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Amp_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_Amp_th as text
%        str2double(get(hObject,'String')) returns contents of tag_Amp_th as a double


% --- Executes during object creation, after setting all properties.
function tag_Amp_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_Amp_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_refresh.
function tag_refresh_Callback(hObject, eventdata, handles)
% hObject    handle to tag_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


A = questdlg('This will overwtite manual annotation. Poceed anyway?');
if ~isequal(A,'Yes')
    return
end
set(handles.tag_refresh,'backgroundColor',[1 0.9 0.4]);
pause(.5);

dat0 =  handles.tag_table.Data;
% Re-sort the table
[~,iis] = sort([dat0{2:end,2}],'ascend');
dat0(2:end,:) = dat0(iis+1,:);

% old.chan_deleted.Amplitude = mean(abs(handles.data.S)>handles.data.ParamPreProc.amp_th).'>0.05;
% old.chan_deleted.SNR = handles.data.SNR(:) < handles.data.ParamPreProc.SNR_th;
% old.chan_deleted.Sig_Corr = mean(handles.data.sig_corr<handles.data.ParamPreProc.sig_corr_th,2)>0.66;
% ii_chan_ko_old = [old.chan_deleted.Amplitude | old.chan_deleted.SNR | old.chan_deleted.Sig_Corr];
% ii_manually_del_now = [dat0{2:end,5}]';
% iidiff = ii_chan_ko_old~=ii_manually_del_now; % these won't be changes because manually modified

% update
handles.data.ParamPreProc.amp_th = str2double(handles.tag_Amp_th.String);
handles.data.ParamPreProc.SNR_th = str2double(handles.tag_SNR_th.String);
handles.data.ParamPreProc.sig_corr_th = str2double(handles.tag_sig_corr_th.String);
handles.data.chan_deleted.Amplitude = mean(abs(handles.data.S)>handles.data.ParamPreProc.amp_th).'>0.05;
handles.data.chan_deleted.SNR = handles.data.SNR(:) < handles.data.ParamPreProc.SNR_th;
handles.data.chan_deleted.Sig_Corr = mean(handles.data.sig_corr<handles.data.ParamPreProc.sig_corr_th,2)>0.66;
% -
handles.data.chan_deleted.ii_chan_ko = [handles.data.chan_deleted.Amplitude | handles.data.chan_deleted.SNR | handles.data.chan_deleted.Sig_Corr];
% handles.data.chan_deleted.ii_chan_ko(iidiff) = ii_manually_del_now(iidiff);
clear iidiff
dat = dat0;
dat(2:end,5) = num2cell(handles.data.chan_deleted.ii_chan_ko);

% =
if isfield(handles.data,'Markers')&&isfield(handles.data,'spikes');
    
    dt = handles.data.Markers.dt - handles.data.spikes(:)*ones(1,size(handles.data.Markers.dt,2));
    rt = handles.data.Markers.rt_Wyatt - handles.data.spikes(:)*ones(1,size(handles.data.Markers.rt_Wyatt,2));
    
    %     ii_chan_check_old = mad(rt,1)*1.48./nanmedian(rt,1)>handles.data.ParamPreProc.cv_thr | mad(dt,1)*1.48./nanmedian(dt,1)>handles.data.ParamPreProc.cv_thr ;
    %     ii_chan_now = [dat0{2:end,6}]';
    %     iidiff = ii_chan_check_old(:) ~= ii_chan_now(:);
    
    handles.data.chan_deleted.ii_chan_check = mad(rt,1)*1.48./nanmedian(rt,1)>str2double(handles.tag_cv_thr.String) | mad(dt,1)*1.48./nanmedian(dt,1)>str2double(handles.tag_cv_thr.String);
else
    handles.data.chan_deleted.ii_chan_check = zeros(size(handles.data.S,2),1);
end

handles.data.ParamPreProc.cv_thr = str2double(handles.tag_cv_thr.String);
dat(2:end,6) = num2cell(handles.data.chan_deleted.ii_chan_check);

handles.tag_table.Data = dat;

set(handles.tag_n_bad_chan,'string',num2str(sum(handles.data.chan_deleted.ii_chan_ko)));
set(handles.tag_n_ok_chan,'string',num2str(sum(~handles.data.chan_deleted.ii_chan_ko)));
set(handles.tag_n_check,'string',num2str(sum(handles.data.chan_deleted.ii_chan_check)));

set(handles.tag_refresh,'backgroundColor',[1 1 1]*.94);
handles.tag_sort_type.Value = 1;
handles.dat = dat;

%% Update
guidata(hObject,handles);


% --- Executes on button press in tag_cv_thr.
function tag_cv_thr_Callback(hObject, eventdata, handles)
% hObject    handle to tag_cv_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function tag_cv_thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_cv_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
