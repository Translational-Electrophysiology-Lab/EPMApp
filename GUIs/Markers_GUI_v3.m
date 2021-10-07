function varargout = Markers_GUI_v3(varargin)
% MARKERS_GUI_V3 MATLAB code for Markers_GUI_v3.fig
%      MARKERS_GUI_V3, by itself, creates a new MARKERS_GUI_V3 or raises the existing
%      singleton*.
%
%      H = MARKERS_GUI_V3 returns the handle to a new MARKERS_GUI_V3 or the handle to
%      the existing singleton*.
%
%      MARKERS_GUI_V3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARKERS_GUI_V3.M with the given input arguments.
%
%      MARKERS_GUI_V3('Property','Value',...) creates a new MARKERS_GUI_V3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Markers_GUI_v3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Markers_GUI_v3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Markers_GUI_v3

% Last Modified by GUIDE v2.5 09-Mar-2018 08:48:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Markers_GUI_v3_OpeningFcn, ...
    'gui_OutputFcn',  @Markers_GUI_v3_OutputFcn, ...
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

% --- Executes just before Markers_GUI_v3 is made visible.
function Markers_GUI_v3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Markers_GUI_v3 (see VARARGIN)

% Choose default command line output for Markers_GUI_v3
handles.output = hObject;
%% ADD DATA
% handles.hObject_main = varargin{1};

vv = fieldnames(varargin{1});
for iv = 1:length(vv)
    handles.(vv{iv}) = varargin{1}.(vv{iv});
end

set(handles.tag_RTmax,'string','450')
set(handles.tag_RTmin,'string','200')
set(handles.tag_DTmax,'string','120')
set(handles.tag_DeltaT,'string','20')
set(handles.tag_LFc,'string','0.5')
set(handles.tag_HFc,'string','25')
set(handles.tag_LFc_AT,'string','0')
set(handles.tag_HFc_AT,'string','100')

if isfield(handles,'MarkersC')
    if isfield(handles.MarkersC,'ParamIn')
        set(handles.tag_RTmax,'string',num2str(handles.MarkersC.ParamIn.max_RT))
        set(handles.tag_RTmin,'string',num2str(handles.MarkersC.ParamIn.min_RT))
        set(handles.tag_DTmax,'string',num2str(handles.MarkersC.ParamIn.DTmax))
        set(handles.tag_LFc,'string',num2str(handles.MarkersC.ParamIn.BW(1)))
        set(handles.tag_HFc,'string',num2str(handles.MarkersC.ParamIn.BW(2)))
        if isfield(handles.MarkersC.ParamIn,'DeltaT')
            set(handles.tag_DeltaT,'string',num2str(handles.MarkersC.ParamIn.DeltaT))
        else
            set(handles.tag_DeltaT,'string','15')
        end
        if isfield(handles.MarkersC.ParamIn,'BW_AT')
            set(handles.tag_LFc_AT,'string',num2str(handles.MarkersC.ParamIn.BW_AT(1)))
            set(handles.tag_HFc_AT,'string',num2str(handles.MarkersC.ParamIn.BW_AT(2)))
        else
            handles.tag_LFc_AT.Visible = 'off';
            handles.tag_HFc_AT.Visible = 'off';
        end
    end
end

set(handles.tag_slider,'sliderStep',[0.050 .2]/size(handles.S,1)*handles.ParamSig.frequency);
% handles.ParamSig = data.ParamSig;
% handles.LOG = data.LOG;
% handles.tag_Log = data.tag_Log;
set([handles.tag_ax_sig,handles.tag_ax_spikes],'xlim',[0 1],'ylim',[0 1],'xtick',[],'ytick',[])

List = handles.ParamSig.Label;
for k = 1:length(List)
    List{k} = ['[',num2str(k),']',List{k}];
end
if isfield(handles,'SNR')
    %     List_new=cell(1,length(List));
    for j=1:length(handles.SNR)
        List{j} = [List{j},'(',num2str(handles.SNR(j),2),'dB)'];
    end
end
if isfield(handles,'sig_corr')
    if ~isempty(handles.sig_corr)
        for j=1:length(List)
            List{j} = [List{j},'/[',num2str(handles.sig_corr(j),2),']'];
        end
    end
end
set(handles.tag_List,'string',List)
clear List*
%            if ~isempty(handles.sig_corr)
%                 List_new{j} = [List{j},'(',num2str(handles.SNR(j),2),'dB)/[',num2str(handles.sig_corr(j),2),']'];
%             else

if ~isfield(handles,'signals_proc_AT');
    handles.tag_show_filter_AT.Visible = 'off';
end
handles.tag_show_filter_RT.Value = 0;
handles.tag_show_filter_RAW.Value = 1;
handles.tag_show_filter_AT.Value = 0;

%%
% data.MarkersC.m_function = 'CorrectFidutialPoints_CathLab.m';
% Update handles structure
guidata(hObject, handles);
% set(handles.tag_List,'string',handles.ParamSig.Label)
% add path for colormaps

% UIWAIT makes Markers_GUI_v3 wait for user response (see UIRESUME)
% uiwait(handles.tag_figure_markers);


% --- Outputs from this function are returned to the command line.
function varargout = Markers_GUI_v3_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% varargout{1} = handles.spikes;
% varargout{1} = getappdata(hObject,'spikes');


% --- Executes on selection change in tag_List.
function tag_List_Callback(hObject, eventdata, handles)
% hObject    handle to tag_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns tag_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_List

tag_Plot_Callback(hObject, [], handles)


% --- Executes during object creation, after setting all properties.
function tag_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_Refresh.
% function tag_Refresh_Callback(hObject, eventdata, handles)
% % hObject    handle to tag_Refresh (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% hold(handles.tag_ax_spikes,'off')
% hold(handles.tag_ax_sig,'off')

% --- Executes on button press in tag_change_polarity_RT.
function tag_change_polarity_RT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_change_polarity_RT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ichan = get(handles.tag_List,'value');
handles.MarkersC.iiTwpos(ichan) = ~handles.MarkersC.iiTwpos(ichan);
% tag_modify_Callback(hObject, [], handles)
%% update
guidata(hObject, handles);

str = cell(1,length(ichan));
for i = 1:length(ichan);
    str{i} = ['IC-',num2str(ichan(i)),': TW polarity from ',num2str(~handles.MarkersC.iiTwpos(ichan(i))),' to ',num2str(handles.MarkersC.iiTwpos(ichan(i)))];
end
msgbox(str);


% --- Executes on button press in tag_modify.
function tag_modify_Callback(hObject, eventdata, handles)
% hObject    handle to tag_modify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


dcm_obj = datacursormode(handles.tag_figure_markers);
set(dcm_obj,'updateFcn',[],'enable','off')
info_struct = getCursorInfo(dcm_obj);
delete(findall(handles.tag_figure_markers,'Type','hggroup','HandleVisibility','off'));
DeltaT = str2double(get(handles.tag_DeltaT,'string')); % ms

if ~isempty(info_struct)
    x = info_struct.Position(1); % [ms]
    xest2 = x-handles.spikes;
    jxest=find(xest2>-DeltaT/1000*handles.ParamSig.frequency);
    jxest = jxest(end);
    %     xest= x-handles.spikes(jxest)+DeltaT; % [ms]
    xest= x-handles.spikes(jxest); % [ms]
else
    xest = nan;
end



ichan = get(handles.tag_List,'value');
if isfield(handles,'SNR')
    snrth_L = str2double(get(handles.tag_SNRth_L,'string'));
    snrth_U = str2double(get(handles.tag_SNRth_U,'string'));
    
    ichan(handles.SNR(ichan)<snrth_L | handles.SNR(ichan)>snrth_U)=[];
    if isempty(ichan)
        helpdlg(['No channels with SNR in [',num2str(snrth_L),'-',num2str(snrth_U),' dB]'])
        return
    end
    set(handles.tag_List,'value',ichan);
end

if get(handles.tag_show_filter_RT,'value')
    signal = handles.signals_proc;
elseif get(handles.tag_show_filter_RAW,'value')
    signal = handles.S;
elseif get(handles.tag_show_filter_AT,'value')
    signal = handles.signals_proc_AT;
end



% Normalization not necessary
% if get(handles.tag_norm,'value')
%     signal = signal./repmat(max(abs(signal)),[size(signal,1) 1]);
% end
%%
% WinMod = round(str2double(get(handles.tag_W_mod,'string'))/1000*handles.ParamSig.frequency); % samples
WinMod = str2double(get(handles.tag_W_mod,'string')); % ms

%%
%% Type of signals
DTmax = str2double(get(handles.tag_DTmax,'string'));
min_RT = str2double(get(handles.tag_RTmin,'string'));
max_RT = str2double(get(handles.tag_RTmax,'string'));

if isfield(handles.MarkersC,'ParamIn')
    ParamIn_C = handles.MarkersC.ParamIn;
else
    ParamIn_C.BW = [0.5000 25];
    ParamIn_C.Twidth_replace_spk = 15;
    ParamIn_C.Thresh_replace_spk = 10;
    ParamIn_C.do_notch = 0;
    ParamIn_C.frequency = handles.ParamSig.frequency;
end
ParamIn_C.DeltaT = DeltaT;
ParamIn_C.DTmax = DTmax;
ParamIn_C.min_RT = min_RT;
ParamIn_C.max_RT = max_RT;

%% Beat to modify
if handles.tag_all_beats.Value
    Nbeats = [1:length(handles.spikes)];
elseif get(handles.tag_ALL_beat_CL,'value')
    % all the heartbeats with similar CL, because if the CL change the
    % morphology of the EGM changes as well
    Nbeats = find(diff(handles.spikes)< handles.spikes(jxest+1)-handles.spikes(jxest)+8 & diff(handles.spikes)> handles.spikes(jxest+1)-handles.spikes(jxest)-8);
    Nbeats = [Nbeats(:)' Nbeats(end)+1]; % need to be a line-vector
    %     T = [1:size(signal,1)];
else
    Nbeats =  jxest;
    
end
% =
if ~get(handles.tag_all_chan,'value')
    ix = round(info_struct.Position(1)*handles.ParamSig.frequency/1000);
    if get(handles.tag_show_filter_RT,'value')
        V = handles.signals_proc;
    elseif get(handles.tag_show_filter_RAW,'value')
        V = handles.S;
    else get(handles.tag_show_filter_AT,'value')
        V = handles.signals_proc_AT;
    end
    if get(handles.tag_dV,'value')
        V=diff(V);
    end
    if get(handles.tag_norm,'value')
        V = V./repmat(max(abs(V)),[size(V,1) 1]);
    end
    
    ichan = (find(V(ix,:)==info_struct.Position(2)));
    clear V
end


Info_Correct.MarkersC = handles.MarkersC;
Info_Correct.xest = xest;
Info_Correct.WinMod = WinMod;
Info_Correct.ichan = ichan;
Info_Correct.Nbeats = Nbeats;
switch handles.tag_modify_type.Value
    case 1
        Info_Correct.Marker_name = 'AT';
    case 2
        if handles.tag_method.Value==2
            Info_Correct.Marker_name = 'RT_Wyatt';
        elseif handles.tag_method.Value==3
            Info_Correct.Marker_name = 'RT_Alt';
        end
    case 3
        Info_Correct.Marker_name = 'ISO';
    case 4
        Info_Correct.Marker_name = 'Tpeak';
    case 5
        Info_Correct.Marker_name = 'Tend';
    case 6
        Info_Correct.Marker_name = 'QRS';
end

do_DTRT_only_fast = 0;
if ~isfield(handles,'MarkersInput');
    if sum(strcmp({'AT','RT_Wyatt','RT_Alt','Tpeak'},Info_Correct.Marker_name))>0
        do_DTRT_only_fast = 1;
        handles.MarkersInput.do_DTRT_only_fast = 1;
    end
else
    if ~isfield(handles.MarkersInput,'do_DTRT_only_fast')
        do_DTRT_only_fast = 1;
        handles.MarkersInput.do_DTRT_only_fast = 1;
    else
        do_DTRT_only_fast = handles.MarkersInput.do_DTRT_only_fast;
    end
end

if do_DTRT_only_fast
    [MarkersC2] = DTRT_mo_gui_v4(signal,handles.spikes,ParamIn_C,Info_Correct);
else
    [MarkersC2] = DTRT_mo_gui_v3(signal,handles.spikes,ParamIn_C,Info_Correct);
end


w1 = xest-WinMod/2;
w2 = xest+WinMod/2;

switch handles.tag_modify_type.Value
    case 1
        A = handles.MarkersC.dt(Nbeats,ichan);
        Anew = MarkersC2.dt;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.dt(Nbeats,ichan) = Anew;
        
        w1 = max(w1,str2double(-handles.tag_DeltaT.String));
        w2 = min(w2,str2double(handles.tag_DTmax.String));
        T = Anew; % for plotting purposes
    case 2
        
        A = handles.MarkersC.rt_down(Nbeats,ichan);
        Anew = MarkersC2.rt_down;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.rt_down(Nbeats,ichan) = Anew;
        
        A = handles.MarkersC.rt_up(Nbeats,ichan);
        Anew = MarkersC2.rt_up;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.rt_up(Nbeats,ichan) = Anew;
        
        if handles.tag_method.Value==3
            A = handles.MarkersC.rt_Alternative(Nbeats,ichan);
            Anew = MarkersC2.rt_Alternative;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.rt_Alternative(Nbeats,ichan) = Anew;
            T = Anew; % for plotting purposes
        elseif handles.tag_method.Value==2
            A = handles.MarkersC.rt_Wyatt(Nbeats,ichan);
            Anew = MarkersC2.rt_Wyatt;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.rt_Wyatt(Nbeats,ichan) = Anew;
            T = Anew; % for plotting purposes
        end
        
        if isfield(handles.MarkersC,'ATw');
            A = handles.MarkersC.ATw(Nbeats,ichan);
            Anew = MarkersC2.ATw;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.ATw(Nbeats,ichan) = Anew;
        end
        
        if isfield(handles.MarkersC,'ATw2');
            A = handles.MarkersC.ATw2(Nbeats,ichan);
            Anew = MarkersC2.ATw2;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.ATw2(Nbeats,ichan) = Anew;
        end
        
        handles.MarkersC.iiTwpos(ichan) = MarkersC2.iiTwpos;
        
        w1 = max(w1,str2double(handles.tag_RTmin.String));
        w2 = min(w2,str2double(handles.tag_RTmax.String));
        
    case 3
        A = handles.MarkersC.isot(Nbeats,ichan);
        Anew = MarkersC2.isot;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.isot(Nbeats,ichan) = Anew;
        
        w1 = max(w1,0);
        w2 = min(w2,120+str2double(handles.tag_DTmax.String));
        T = Anew; % for plotting purposes
    case 4
        A = handles.MarkersC.tTpeak(Nbeats,ichan);
        Anew = MarkersC2.tTpeak;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.tTpeak(Nbeats,ichan) = Anew;
        w1 = max(w1,str2double(handles.tag_RTmin.String));
        w2 = min(w2,str2double(handles.tag_RTmax.String));
        T = Anew; % for plotting purposes
        
        A = handles.MarkersC.Tpeak_amp(Nbeats,ichan);
        Anew = MarkersC2.Tpeak_amp;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.Tpeak_amp(Nbeats,ichan) = Anew;
        
    case 5
        A = handles.MarkersC.tTend(Nbeats,ichan);
        Anew = MarkersC2.tTend;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.tTend(Nbeats,ichan) = Anew;
        w1 = max(w1,str2double(handles.tag_RTmin.String));
        w2 = min(w2,str2double(handles.tag_RTmax.String));
        T = Anew;
        
        if isfield(handles.MarkersC,'tTend_defl');
            A = handles.MarkersC.tTend_defl(Nbeats,ichan);
            Anew = MarkersC2.tTend_defl;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.tTend_defl(Nbeats,ichan) = Anew;
        end
        
        A = handles.MarkersC.rt_Alternative(Nbeats,ichan);
        Anew = MarkersC2.rt_Alternative;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.rt_Alternative(Nbeats,ichan) = Anew;
        
        A = handles.MarkersC.rt_down(Nbeats,ichan);
        Anew = MarkersC2.rt_down;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.rt_down(Nbeats,ichan) = Anew;
        
        A = handles.MarkersC.rt_up(Nbeats,ichan);
        Anew = MarkersC2.rt_up;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.rt_up(Nbeats,ichan) = Anew;
        
        
        w1 = max(w1,str2double(handles.tag_RTmin.String));
        w2 = min(w2,str2double(handles.tag_RTmax.String));
        
    case 6
        A = handles.MarkersC.tQRSoff(Nbeats,ichan);
        Anew = MarkersC2.tQRSoff;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.tQRSoff(Nbeats,ichan) = Anew;
        clear A Anew
        
        A = handles.MarkersC.tQRSon(Nbeats,ichan);
        Anew = MarkersC2.tQRSon;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.tQRSon(Nbeats,ichan) = Anew;
        clear A Anew
        
        A = handles.MarkersC.tQw(Nbeats,ichan);
        Anew = MarkersC2.tQw;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.tQw(Nbeats,ichan) = Anew;
        clear A Anew
        
        A = handles.MarkersC.tRw(Nbeats,ichan);
        Anew = MarkersC2.tRw;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.tRw(Nbeats,ichan) = Anew;
        T = Anew; % for plotting purposes
        clear A Anew
        
        A = handles.MarkersC.tSw(Nbeats,ichan);
        Anew = MarkersC2.tSw;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.tSw(Nbeats,ichan) = Anew;
        clear A Anew
        
        A = handles.MarkersC.QRS_amp(Nbeats,ichan);
        Anew = MarkersC2.QRS_amp;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.QRS_amp(Nbeats,ichan) = Anew;
        clear A Anew
        
        
        A = handles.MarkersC.QRS_area(Nbeats,ichan);
        Anew = MarkersC2.QRS_area;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.QRS_area(Nbeats,ichan) = Anew;
        clear A Anew
        
        A = handles.MarkersC.QRSoff(Nbeats,ichan);
        Anew = MarkersC2.QRSoff;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.QRSoff(Nbeats,ichan) = Anew;
        clear A Anew
        
        
        A = handles.MarkersC.QRSon(Nbeats,ichan);
        Anew = MarkersC2.QRSon;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.QRSon(Nbeats,ichan) = Anew;
        clear A Anew
        
        A = handles.MarkersC.QRSw(Nbeats,ichan);
        Anew = MarkersC2.QRSw;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.QRSw(Nbeats,ichan) = Anew;
        clear A Anew
        
        A = handles.MarkersC.QRSw_fw90(Nbeats,ichan);
        Anew = MarkersC2.QRSw_fw90;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
        handles.MarkersC.QRSw_fw90(Nbeats,ichan) = Anew;
        clear A Anew
        
        w1 = max(w1,str2double(-handles.tag_DeltaT.String));
        w2 = min(w2,str2double(handles.tag_DTmax.String));
end
clear A*
tag_Plot_Callback(hObject, [], handles);

% if handles.tag_all_beats.Value
% %     xs = find(xest>handles.spikes);xs = handles.spikes-handles.spikes(is);
%     for i = 1:length(Nbeats)
% %     plot(handles.tag_ax_sig,xs(i)+w1*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
% %     plot(handles.tag_ax_sig,xs(i)+w2*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
%     plot(handles.tag_ax_sig,handles.spikes(Nbeats(i))+w1*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
%     plot(handles.tag_ax_sig,handles.spikes(Nbeats(i))+w2*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
%
%     end
% else
%     plot(handles.tag_ax_sig,handles.spikes(Nbeats)+w1*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
%     plot(handles.tag_ax_sig,handles.spikes(Nbeats)+w2*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
% end
%

% normalize, just for visualization
if get(handles.tag_norm,'value')
    signal = signal./repmat(max(abs(signal)),[size(signal,1) 1]);
end
T = round(nanmedian(T,2)/1000*handles.ParamSig.frequency);
Am = zeros(size(Nbeats));
Am(~isnan(T)) = nanmedian(signal(T(~isnan(T)),ichan),2);

% for i = 1:length(Nbeats)
%     plot(handles.tag_ax_sig,handles.spikes(Nbeats(i))+w1*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+Am(i),'--r')
%     plot(handles.tag_ax_sig,handles.spikes(Nbeats(i))+w2*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+Am(i),'--r')
% end
plot(handles.tag_ax_sig,handles.spikes(Nbeats)'+w1*[1;1],[-1;1]*range(get(handles.tag_ax_sig,'ylim'))/3+Am,'--r')
plot(handles.tag_ax_sig,handles.spikes(Nbeats)'+w2*[1;1],[-1;1]*range(get(handles.tag_ax_sig,'ylim'))/3+Am,'--r')

%%
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in tag_Plot.
function tag_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.tag_hold.Value
    hold(handles.tag_ax_sig,'on');
    hold(handles.tag_ax_spikes,'on');
else
    hold(handles.tag_ax_sig,'off');
    hold(handles.tag_ax_spikes,'off');
end

if isfield(handles,'MarkersC')
    v= fieldnames(handles.MarkersC);
    for i=1:length(v)
        eval([v{i},'=handles.MarkersC.',(v{i}),';']);
    end
else
    v= fieldnames(handles.Markers);
    for i=1:length(v)
        eval([v{i},'=handles.Markers.',(v{i}),';']);
    end
end

if ~exist('ParamIn','var')
    ParamIn = handles.ParamSig;
end

if ~isfield(ParamIn,'frequency')
    ParamIn.frequency = handles.ParamSig.frequency;
end
if length(ParamIn)>1
    ParamIn = ParamIn(1);
end
ari_Wyatt = rt_Wyatt-dt;
ari_Alternative = rt_Alternative-dt;

if get(handles.tag_method,'value')==1
    set(handles.tag_method,'value',2)
end
iiv = get(handles.tag_List,'value');
% show only channels within SNRth_L and SNR_U
if isfield(handles,'SNR')
    snrth_L = str2double(get(handles.tag_SNRth_L,'string'));
    snrth_U = str2double(get(handles.tag_SNRth_U,'string'));
    
    iiv(handles.SNR(iiv)<snrth_L | handles.SNR(iiv)>snrth_U)=[];
    if isempty(iiv)
        helpdlg(['No channels with SNR in [',num2str(snrth_L),'-',num2str(snrth_U),' dB]'])
        return
    end
    set(handles.tag_List,'value',iiv);
end


% Filtered/not filtered
if get(handles.tag_show_filter_RAW,'value')
    signals=handles.S;
elseif get(handles.tag_show_filter_RT,'value')
    signals=handles.signals_proc;
elseif get(handles.tag_show_filter_AT,'value')
    signals=handles.signals_proc_AT;
else
    h = helpdlg('Please select one of the filter mode: RAW, RT, AT');
    return
end

% Unipolar/Bipolar
if get(handles.tag_plot_bipolar,'value')
    addpath .\GUI_egm_mFiles\Geo_Chann
    if isempty(handles.ParamSig.geoname)|isequal(handles.ParamSig.geoname,'Unknown')
        geo = inputdlg('Please enter name of sock geometry','Geometry',1,{'mo_sock2'})
        handles.ParamSig.geoname = geo{1};
        geo = load(['ALLgeoDATA_',geo{1}]);
    else
        geo = load(['ALLgeoDATA_',handles.ParamSig.geoname]);
    end
    s2 = nan(size(signals));
    for j = 1:size(signals,2)
        if ~isnan( geo.channels_closes.ref_bipol(j))&geo.channels_closes.ref_bipol(j)<=size(signals,2)
            s2(:,j) = signals(:,j) - signals(:,geo.channels_closes.ref_bipol(j));
        end
    end
    set(handles.tag_dV,'value',0)
    signals = s2;
end

% dV/dt
if get(handles.tag_dV,'value')
    signals=diff(signals);
end

if get(handles.tag_norm,'value')
    signals = signals./repmat(max(abs(signals)),[size(signals,1) 1]);
end

% xx = xlim(handles.tag_ax_sig);
% aa = plot(handles.tag_ax_sig,[1:size(signals,1)]/ParamIn.frequency*1000,signals(:,iiv));

% for Matlab2015
xx = xlim(handles.tag_ax_sig);
aa = plot(handles.tag_ax_sig,[1:size(signals,1)]/ParamIn.frequency*1000,signals(:,iiv));
hold(handles.tag_ax_sig,'on')
zoom(handles.tag_ax_spikes,'reset')
zoom(handles.tag_ax_sig,'reset')

% Set color manualy
% if length(aa)>1
%     cmap = jet(length(aa));
%     icmap = find(cmap(:,1)==1&cmap(:,2)==1&cmap(:,3)==0);
%     if ~isempty(icmap)
%         cmap(icmap,:)=[0 0 0]; % yellow->black
%     end
%     %     cmap(cmap(:,1)==1&cmap(:,2)==1&cmap(:,3)==0,:)=[0 0 0]; % yellow->black
% else
%     cmap = [0 0 1];
% end
% for i =1:length(aa)
%     set(aa(i),'color',cmap(i,:))
% end
% % ==


yy =get(handles.tag_ax_sig,'ylim');
yym = max(.99*max(yy));yym2 = max(1.01*min(yy));
plot(handles.tag_ax_sig,handles.spikes,yym,'vk');
plot(handles.tag_ax_sig,handles.spikes,yym2,'^k');
if length(handles.spikes)<500
    for j=1:length(handles.spikes)
        plot(handles.tag_ax_sig,[1 1]*handles.spikes(j),[yym yym2],':k')
    end
end
Leg = cell(1,length(iiv));
if ~exist('tRw','var')
    tRw = nan(size(dt));
end
for j = 1:length(iiv)
    ichan=iiv(j);
    dtp = dt(:,ichan);dtp(isnan(dtp))=[];
    tRwp = tRw(:,ichan);tRwp(isnan(tRwp))=[];
    
    tTp = tTpeak(:,ichan);tTp(isnan(tTp))=[];
    if get(handles.tag_method,'value')==2
        rtp = rt_Wyatt(:,ichan);rtp(isnan(rtp))=[];
    elseif get(handles.tag_method,'value')==3
        rtp = rt_Alternative(:,ichan);rtp(isnan(rtp))=[];
    end
    isotp = isot(:,ichan);isotp(isnan(isotp))=[];
    tTe = tTend(:,ichan);tTe(isnan(tTe))=[];
    
    if ~isempty(dtp)
        dtp(dtp<1/ParamIn.frequency*1000 | dtp>size(signals,1)/ParamIn.frequency*1000)=[];
        ii = round(dtp*ParamIn.frequency/1000);
        adt_sig(j) = plot(handles.tag_ax_sig,dtp,signals(ii,ichan),'marker','o','color',get(aa(j),'color'),'linestyle','none','markerfacecolor',get(aa(j),'color'),'markeredgecolor','k');
    end
    if ~isempty(rtp)
        rtp(rtp>size(signals,1)/ParamIn.frequency*1000)=[];
        ii = round(rtp*ParamIn.frequency/1000);
        art_sig(j) = plot(handles.tag_ax_sig,rtp,signals(ii,ichan),'v','markersize',8,'color',get(aa(j),'color'),'linewidth',2,'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k');
    end
    if ~isempty(isotp)
        ii = round(isotp*ParamIn.frequency/1000);
        aiso_sig(j) = plot(handles.tag_ax_sig,isotp,signals(ii,ichan),'square','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k');
    end
    if ~isempty(tTp)
        ii = round(tTp*ParamIn.frequency/1000);
        atp_sig(j) = plot(handles.tag_ax_sig,tTp,signals(ii,ichan),'Diamond','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k');
    end
    if ~isempty(tTe)
        ii = round(tTe*ParamIn.frequency/1000);tTe(ii<1)=[];ii(ii<1)=[];
        tTe(ii>size(signals,1))=[];ii(ii>size(signals,1))=[];
        if ~isempty(ii)
            atend_sig(j) = plot(handles.tag_ax_sig,tTe,signals(ii,ichan),'^','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k');
        end
    end
    if ~isempty(tRwp)
        ii = round(tRwp*ParamIn.frequency/1000);
        atRw_sig(j) = plot(handles.tag_ax_sig,tRwp,signals(ii,ichan),'Diamond','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k');
    end
    
    clear ii
    Leg(j) =  {[handles.ParamSig.Label{ichan},'(',num2str(handles.SNR(ichan),2),'dB)']};
end
ll = legend(aa,Leg,'AutoUpdate','off');
if ~get(handles.tag_legend,'value')
    set(ll,'visible','off')
end
% --
for j=1:length(iiv)
    adt(j)=plot(handles.tag_ax_spikes,dt(:,iiv(j)),dt(:,iiv(j))-handles.spikes);
    hold(handles.tag_ax_spikes,'on')
    if get(handles.tag_method,'value')==2
        aari(j)=plot(handles.tag_ax_spikes,rt_Wyatt(:,iiv(j)),ari_Wyatt(:,iiv(j)));
        art(j)=plot(handles.tag_ax_spikes,rt_Wyatt(:,iiv(j)),rt_Wyatt(:,iiv(j))-handles.spikes);
    elseif get(handles.tag_method,'value')==3
        aari(j)=plot(handles.tag_ax_spikes,rt_Alternative(:,iiv(j)),ari_Alternative(:,iiv(j)));
        art(j)=plot(handles.tag_ax_spikes,rt_Alternative(:,iiv(j)),rt_Alternative(:,iiv(j))-handles.spikes);
    end
    
    aiso(j)=plot(handles.tag_ax_spikes,isot(:,iiv(j)),isot(:,iiv(j))-handles.spikes);
    atp(j)=plot(handles.tag_ax_spikes,tTpeak(:,iiv(j)),tTpeak(:,iiv(j))-handles.spikes);
    ate(j)=plot(handles.tag_ax_spikes,tTend(:,iiv(j)),tTend(:,iiv(j))-handles.spikes);
    atR(j)=plot(handles.tag_ax_spikes,tRw(:,iiv(j)),tRw(:,iiv(j))-handles.spikes);
    if size(tRw,1) > 1
        aRRI(j)=plot(handles.tag_ax_spikes,tRw(2:end,iiv(j)),diff(tRw(:,iiv(j))));
        set(aRRI(j),'marker','o','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k','linestyle','--')
    end
    
    set(adt(j),'marker','o','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k','linestyle',':')
    set(aari(j),'marker','*','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k','linestyle',':')
    set(art(j),'marker','V','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k','linestyle',':')
    set(aiso(j),'marker','square','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k','linestyle',':')
    set(atp(j),'marker','diamond','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k','linestyle',':')
    set(ate(j),'marker','^','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k','linestyle',':')
    set(atR(j),'marker','square','color',get(aa(j),'color'),'markerfacecolor',get(aa(j),'color'),'markeredgecolor','k','linestyle',':')
end


if ~get(handles.tag_show_DT,'value')
    if exist('adt_sig','var')
        set(adt_sig(isgraphics(adt_sig)),'visible','off')
    end
    set(adt,'visible','off')
end
if ~get(handles.tag_show_RT,'value')
    if exist('art_sig','var')
        set(art_sig(isgraphics(art_sig)),'visible','off')
    end
    set(art,'visible','off')
end
if ~get(handles.tag_show_ARI,'value')
    set(aari(isgraphics(aari)),'visible','off')
end
if ~get(handles.tag_show_Tiso,'value')
    if exist('aiso_sig','var')
        set(aiso_sig(isgraphics(aiso_sig)),'visible','off')
    end
    set(aiso,'visible','off')
end
if ~get(handles.tag_show_Tpeak,'value')
    if exist('atp_sig','var')
        set(atp_sig(isgraphics(atp_sig)),'visible','off')
    end
    set(atp,'visible','off')
end
if ~get(handles.tag_show_Tend,'value')
    if exist('atend_sig','var')
        set(atend_sig(isgraphics(atend_sig)),'visible','off')
    end
    set(ate,'visible','off')
end
if ~get(handles.tag_show_Rw,'value')
    if exist('atRw_sig','var')
        set(atRw_sig(isgraphics(atRw_sig)),'visible','off')
    end
    set(atR,'visible','off')
end

if get(handles.tag_show_RRI,'value')
    if exist('atRw_sig','var')
        set(atRw_sig,'visible','on')
    end
else
    if exist('aRRI','var')
        set(aRRI,'visible','off')
    end
end

if size(tRw,1)>1
    ll2=legend([adt(1),aiso(1),aari(1),art(1),atp(1),ate(1),atR(1),aRRI(1)],'DT','ISO','ARI','RT','Tp','Tend','Rw','RRI');
else
    ll2=legend([adt(1),aiso(1),aari(1),art(1),atp(1),ate(1),atR(1)],'DT','ISO','ARI','RT','Tp','Tend','Rw');
end
linkaxes([handles.tag_ax_spikes,handles.tag_ax_sig],'x')

if ~isequal(xx,[0 1])&xx(2)>200
    set([handles.tag_ax_sig],'xlim',xx,'ylim',yy)
    set([handles.tag_ax_spikes],'xlim',xx)
end

set(handles.tag_ax_sig,'xticklabel',[])
set(handles.tag_ax_spikes,'ygrid','on')
xlabel(handles.tag_ax_spikes,'(ms)')
handles.yym = yym;
ylabel(handles.tag_ax_spikes,'[ms]')
% update
guidata(hObject, handles);

% keyboard
%
dcm_obj = datacursormode(gcf);
% set(dcm_obj,'UpdateFcn',@myfunctioncursor_markers)


% --- Executes on button press in tag_done.
function tag_done_Callback(hObject, eventdata, handles)
% hObject    handle to tag_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = guidata(handles.hObject_main);

if isfield(data,'Resti')
    answ = questdlg('Eliminate Restitution based on old markers?','Restitution','Yes','No','Yes');
    
    if isequal(answ,'Yes')
        data=rmfield(data,'Resti');
    end
end

data.MarkersC = handles.MarkersC;
if isfield(handles,'Markers')
    data.Markers = handles.Markers;
end
data.SNR = handles.SNR;
if isfield(handles,'sig_corr')
    data.sig_corr = handles.sig_corr;
end

if isfield('handles','sig_corr_beats')
    data.sig_corr_beats = handles.sig_corr_beats;
end

data.signals_proc = handles.signals_proc;

if isfield(handles,'signals_proc_AT');
    data.signals_proc_AT = handles.signals_proc_AT;
end
data.ParamSig = handles.ParamSig;
% update
guidata(handles.hObject_main,data);
close(handles.tag_figure_markers)



function tag_SNRth_L_Callback(hObject, eventdata, handles)
% hObject    handle to tag_SNRth_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_SNRth_L as text
%        str2double(get(hObject,'String')) returns contents of tag_SNRth_L as a double


% --- Executes during object creation, after setting all properties.
function tag_SNRth_L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_SNRth_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function tag_tools_Callback(hObject, eventdata, handles)
% hObject    handle to tag_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in tag_method.
function tag_method_Callback(hObject, eventdata, handles)
% hObject    handle to tag_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_method


% --- Executes during object creation, after setting all properties.
function tag_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tag_modify_type.
function tag_modify_type_Callback(hObject, eventdata, handles)
% hObject    handle to tag_modify_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_modify_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_modify_type


% --- Executes during object creation, after setting all properties.
function tag_modify_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_modify_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_all_beats.
function tag_all_beats_Callback(hObject, eventdata, handles)
% hObject    handle to tag_all_beats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_all_beats



function tag_DTmax_Callback(hObject, eventdata, handles)
% hObject    handle to tag_DTmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_DTmax as text
%        str2double(get(hObject,'String')) returns contents of tag_DTmax as a double


% --- Executes during object creation, after setting all properties.
function tag_DTmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_DTmax (see GCBO)
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



function tag_distTospike_Callback(hObject, eventdata, handles)
% hObject    handle to tag_distTospike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_distTospike as text
%        str2double(get(hObject,'String')) returns contents of tag_distTospike as a double


% --- Executes during object creation, after setting all properties.
function tag_distTospike_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_distTospike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in tag_find_markers.
function tag_find_markers_Callback(hObject, eventdata, handles)
% hObject    handle to tag_find_markers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% clear variables which are gona be overrwritten
set(handles.tag_find_markers,'backgroundcolor',[.6 .6 .6],'foregroundcolor','r','FontWeight','bold')


vover = {'signals_raw','signals_proc','S','Markers','MarkersC','SNR'};
for i=1:length(vover)
    if isfield(handles,vover{i})
        handles = rmfield(handles,vover{i});
    end
end

hold(handles.tag_ax_sig,'off')
plot(handles.tag_ax_sig,0,0,'visible','off')
text(.5,.5,'Localizing new markers ...')
hold(handles.tag_ax_spikes,'off')
plot(handles.tag_ax_spikes,0,0,'visible','off')
set([handles.tag_ax_sig,handles.tag_ax_spikes],'xlim',[0 1],'ylim',[0 1],'xtick',[],'ytick',[])
set(handles.tag_List,'value',[1:length(get(handles.tag_List,'string'))])
%%
% hinput = Markers_input(handles.tag_figure_markers);
hinput = Markers_input_v2(handles.tag_figure_markers);

% update
handles = guidata(hObject);
handles.MarkersInput.frequency = handles.ParamSig.frequency; % sampling frequency
if ~handles.MarkersInput.ok
    set(handles.tag_find_markers,'backgroundcolor',[1 1 1]*.941,'foregroundcolor','k','FontWeight','normal')
    return
end
%
spikes = handles.spikes;
spikes_samp = round(spikes/1000*handles.ParamSig.frequency);
if isequal(handles.type_analysis,'Pacing')
    hwb = waitbar(1/4,'Preprocessing [Pacing]...');
else
    hwb = waitbar(1/4,'Preprocessing [Sinus Rhythm]...');
end
% delete spikes
if ~isequal(handles.type_analysis,'Sinus Rhythm')
    [signals_nosp,~] = replace_spikes_mo(handles.S,spikes_samp,round(handles.MarkersInput.Twidth_replace_spk/1000*handles.ParamSig.frequency),handles.MarkersInput.Thresh_replace_spk);
else
    signals_nosp = handles.S;
end


waitbar(2/4,hwb,'... filtering');
%% =======================================
%% Markers (AT in non-filtered - RT in filtered)
ParamIn = handles.MarkersInput;
if isfield(ParamIn,'BW_AT')
    if isfield(handles,'MarkersC');
        %         set(handles.tag_LFc_AT,'string',num2str(handles.MarkersC.ParamIn.BW_AT(1)))
        %         set(handles.tag_HFc_AT,'string',num2str(handles.MarkersC.ParamIn.BW_AT(2)))
        set(handles.tag_LFc_AT,'string',num2str(ParamIn.BW_AT(1)))
        set(handles.tag_HFc_AT,'string',num2str(ParamIn.BW_AT(2)))
    end
    handles.tag_show_filter_AT.Visible = 'on';
else
    handles.tag_LFc_AT.Visible = 'off';
    handles.tag_HFc_AT.Visible = 'off';
    handles.tag_show_filter_AT.Visible = 'off';
end

signals_proc = butterworthfilter(signals_nosp,handles.ParamSig.frequency,ParamIn.BW);
signals_proc_AT = butterworthfilter(signals_nosp,handles.ParamSig.frequency,ParamIn.BW_AT);
if ParamIn.do_notch;
    wo = 50/(ParamIn.frequency/2);
    bw = wo/35;
    [b,a] = iirnotch(wo,bw);
    signals_proc_AT = filtfilt(b,a,signals_proc_AT);
end
% ParamIn.DeltaT = 15;

if handles.MarkersInput.do_DTRT_only_fast
    Markers_TW = DTRT_mo_gui_v4(signals_proc,spikes,ParamIn);
    Markers_QRS = DTRT_mo_gui_v4(signals_proc_AT,spikes,ParamIn);
else
    Markers_TW = DTRT_mo_gui_v3(signals_proc,spikes,ParamIn);
    Markers_QRS = DTRT_mo_gui_v3(signals_proc_AT,spikes,ParamIn);
end

params_from_TW = {'tTpeak','Tpeak_amp','tPwave_int','tPw','tTend','tTend_defl','rt_Wyatt','rt_Alternative','iiTwpos','rt_up','isot','rt_down','ATw','ATw2'};
params_from_QRS = {'tSw','tRw','tQw','tQRSon','tQRSoff','QRSw_fw90','QRSw','dt','tdtM','tdtm','RatioI','ParamIn','m_function','Legend',...,
    'Rw_amp','Sw_amp','Qw_amp','QRS_amp','QRS_area'};

for i = 1:length(params_from_QRS)
    Markers.(params_from_QRS{i}) = Markers_QRS.(params_from_QRS{i});
end
for i = 1:length(params_from_TW)
    Markers.(params_from_TW{i}) = Markers_TW.(params_from_TW{i});
end

%% =======================================
[~,SNR]=snr_mo_fast_few_beats(signals_nosp,handles.spikes,handles.ParamSig.frequency,handles.SNR_Bsig,handles.SNR_Bnoise);
if length(spikes)>1
    [~,sig_corr] = sig_corr_mo_all_beats(signals_proc_AT,spikes_samp);
else
    sig_corr = nan(size(signals_proc,2),1);
end
%%
List = get(handles.tag_List,'string');
List_new=cell(1,length(List));
for j=1:length(List)
    aa = List{j};ii = find(aa=='(');
    if ~isempty(ii)
        aa(ii(1):end)=[];
    end
    if ~isempty(sig_corr)
        List_new{j} = [aa,'(',num2str(SNR(j),2),'dB)/[',num2str(sig_corr(j),2),']'];
    else
        List_new{j} = [aa,'(',num2str(SNR(j),2),'dB)'];
    end
end
set(handles.tag_List,'string',List_new)
clear List*


handles.Markers = Markers;
handles.MarkersC = Markers;
handles.signals_proc = signals_proc;
handles.signals_proc_AT = signals_proc_AT;
handles.SNR = SNR;
handles.sig_corr = sig_corr;
% handles.sig_corr_beats = sig_corr_beats;
% handles.sig_corr100 = sig_corr100;
% handles.sig_corr100_beats = sig_corr100_beats;
%%
set(handles.tag_RTmax,'string',num2str(handles.MarkersInput.max_RT))
set(handles.tag_RTmin,'string',num2str(handles.MarkersInput.min_RT))
set(handles.tag_DTmax,'string',num2str(handles.MarkersInput.DTmax))
set(handles.tag_LFc,'string',num2str(handles.MarkersInput.BW(1)))
set(handles.tag_HFc,'string',num2str(handles.MarkersInput.BW(2)))
set(handles.tag_LFc_AT,'string',num2str(handles.MarkersInput.BW_AT(1)))
set(handles.tag_HFc_AT,'string',num2str(handles.MarkersInput.BW_AT(2)))
set(handles.tag_DeltaT,'string',num2str(handles.MarkersInput.DeltaT))

%% update
guidata(hObject, handles);
set(handles.tag_find_markers,'backgroundcolor',[1 1 1]*.941,'foregroundcolor','k','FontWeight','normal')

close(hwb)


function tag_LFc_AT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_LFc_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_LFc_AT as text
%        str2double(get(hObject,'String')) returns contents of tag_LFc_AT as a double


% --- Executes during object creation, after setting all properties.
function tag_LFc_AT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_LFc_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_LFc_Callback(hObject, eventdata, handles)
% hObject    handle to tag_LFc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_LFc as text
%        str2double(get(hObject,'String')) returns contents of tag_LFc as a double


% --- Executes during object creation, after setting all properties.
function tag_LFc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_LFc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_Autocorrection.
function tag_Autocorrection_Callback(hObject, eventdata, handles)

%%
DeltaT = str2double(get(handles.tag_DeltaT,'string')); % ms
ichan = get(handles.tag_List,'value');
if isfield(handles,'SNR')
    snrth_L = str2double(get(handles.tag_SNRth_L,'string'));
    snrth_U = str2double(get(handles.tag_SNRth_U,'string'));
    
    ichan(handles.SNR(ichan)<snrth_L | handles.SNR(ichan)>snrth_U)=[];
    if isempty(ichan)
        helpdlg(['No channels with SNR in [',num2str(snrth_L),'-',num2str(snrth_U),' dB]'])
        return
    end
    set(handles.tag_List,'value',ichan);
end

if get(handles.tag_show_filter_RT,'value')
    signal = handles.signals_proc;
elseif get(handles.tag_show_filter_RAW,'value')
    signal = handles.S;
elseif get(handles.tag_show_filter_AT,'value')
    signal = handles.signals_proc_AT;
end



WinMod = str2double(get(handles.tag_W_mod,'string')); % ms

%%
%% Type of signals
DTmax = str2double(get(handles.tag_DTmax,'string'));
min_RT = str2double(get(handles.tag_RTmin,'string'));
max_RT = str2double(get(handles.tag_RTmax,'string'));

if isfield(handles.MarkersC,'ParamIn')
    ParamIn_C = handles.MarkersC.ParamIn;
else
    ParamIn_C.BW = [0.5000 25];
    ParamIn_C.Twidth_replace_spk = 15;
    ParamIn_C.Thresh_replace_spk = 10;
    ParamIn_C.do_notch = 0;
    ParamIn_C.frequency = handles.ParamSig.frequency;
end
ParamIn_C.DeltaT = DeltaT;
ParamIn_C.DTmax = DTmax;
ParamIn_C.min_RT = min_RT;
ParamIn_C.max_RT = max_RT;

Nbeats = [1:length(handles.spikes)];

Info_Correct.MarkersC = handles.MarkersC;
% Info_Correct.xest = xest;
Info_Correct.WinMod = WinMod;
% Info_Correct.ichan = ichan;
Info_Correct.Nbeats = Nbeats;
switch handles.tag_modify_type.Value
    case 1
        Info_Correct.Marker_name = 'AT';
        xest_all = nanmedian(handles.MarkersC.dt-handles.spikes,1);
    case 2
        if handles.tag_method.Value==2
            Info_Correct.Marker_name = 'RT_Wyatt';
            xest_all = nanmedian(handles.MarkersC.rt_Wyatt-handles.spikes,1);
        elseif handles.tag_method.Value==3
            Info_Correct.Marker_name = 'RT_Alt';
            xest_all = nanmedian(handles.MarkersC.rt_Alternative-handles.spikes,1);
        end
    case 3
        Info_Correct.Marker_name = 'ISO';
        xest_all = nanmedian(handles.MarkersC.isot-handles.spikes,1);
    case 4
        Info_Correct.Marker_name = 'Tpeak';
        xest_all = nanmedian(handles.MarkersC.tTpeak-handles.spikes,1);
    case 5
        Info_Correct.Marker_name = 'Tend';
        xest_all = nanmedian(handles.MarkersC.tTend-handles.spikes,1);
    case 6
        Info_Correct.Marker_name = 'QRS';
end

ichan_tot = ichan;
clear ichan
for j = 1:length(ichan_tot)
    Info_Correct.xest = xest_all(ichan_tot(j));
    Info_Correct.ichan = ichan_tot(j);
    ichan = ichan_tot(j);
    Info_Correct.MarkersC = handles.MarkersC;
    
    
%     if handles.MarkersInput.do_DTRT_only_fast
        
        [MarkersC2] = DTRT_mo_gui_v4(signal,handles.spikes,ParamIn_C,Info_Correct);
%     else
%         [MarkersC2] = DTRT_mo_gui_v3(signal,handles.spikes,ParamIn_C,Info_Correct);
%     end
%     
    
    w1 = Info_Correct.xest-WinMod/2;
    w2 = Info_Correct.xest+WinMod/2;
    
    switch handles.tag_modify_type.Value
        case 1
            A = handles.MarkersC.dt(Nbeats,ichan);
            Anew = MarkersC2.dt;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.dt(Nbeats,ichan) = Anew;
            
            w1 = max(w1,str2double(-handles.tag_DeltaT.String));
            w2 = min(w2,str2double(handles.tag_DTmax.String));
            T = Anew; % for plotting purposes
        case 2
            
            A = handles.MarkersC.rt_down(Nbeats,ichan);
            Anew = MarkersC2.rt_down;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.rt_down(Nbeats,ichan) = Anew;
            
            A = handles.MarkersC.rt_up(Nbeats,ichan);
            Anew = MarkersC2.rt_up;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.rt_up(Nbeats,ichan) = Anew;
            
            if handles.tag_method.Value==3
                A = handles.MarkersC.rt_Alternative(Nbeats,ichan);
                Anew = MarkersC2.rt_Alternative;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
                handles.MarkersC.rt_Alternative(Nbeats,ichan) = Anew;
                T = Anew; % for plotting purposes
            elseif handles.tag_method.Value==2
                A = handles.MarkersC.rt_Wyatt(Nbeats,ichan);
                Anew = MarkersC2.rt_Wyatt;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
                handles.MarkersC.rt_Wyatt(Nbeats,ichan) = Anew;
                T = Anew; % for plotting purposes
            end
            
            if isfield(handles.MarkersC,'ATw');
                A = handles.MarkersC.ATw(Nbeats,ichan);
                Anew = MarkersC2.ATw;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
                handles.MarkersC.ATw(Nbeats,ichan) = Anew;
            end
            
            if isfield(handles.MarkersC,'ATw2');
                A = handles.MarkersC.ATw2(Nbeats,ichan);
                Anew = MarkersC2.ATw2;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
                handles.MarkersC.ATw2(Nbeats,ichan) = Anew;
            end
            
            handles.MarkersC.iiTwpos(ichan) = MarkersC2.iiTwpos;
            
            w1 = max(w1,str2double(handles.tag_RTmin.String));
            w2 = min(w2,str2double(handles.tag_RTmax.String));
            
        case 3
            A = handles.MarkersC.isot(Nbeats,ichan);
            Anew = MarkersC2.isot;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.isot(Nbeats,ichan) = Anew;
            
            w1 = max(w1,0);
            w2 = min(w2,120+str2double(handles.tag_DTmax.String));
            T = Anew; % for plotting purposes
        case 4
            A = handles.MarkersC.tTpeak(Nbeats,ichan);
            Anew = MarkersC2.tTpeak;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.tTpeak(Nbeats,ichan) = Anew;
            w1 = max(w1,str2double(handles.tag_RTmin.String));
            w2 = min(w2,str2double(handles.tag_RTmax.String));
            T = Anew; % for plotting purposes
            
            A = handles.MarkersC.Tpeak_amp(Nbeats,ichan);
            Anew = MarkersC2.Tpeak_amp;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.Tpeak_amp(Nbeats,ichan) = Anew;
            
        case 5
            A = handles.MarkersC.tTend(Nbeats,ichan);
            Anew = MarkersC2.tTend;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.tTend(Nbeats,ichan) = Anew;
            w1 = max(w1,str2double(handles.tag_RTmin.String));
            w2 = min(w2,str2double(handles.tag_RTmax.String));
            T = Anew;
            
            if isfield(handles.MarkersC,'tTend_defl');
                A = handles.MarkersC.tTend_defl(Nbeats,ichan);
                Anew = MarkersC2.tTend_defl;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
                handles.MarkersC.tTend_defl(Nbeats,ichan) = Anew;
            end
            
            A = handles.MarkersC.rt_Alternative(Nbeats,ichan);
            Anew = MarkersC2.rt_Alternative;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.rt_Alternative(Nbeats,ichan) = Anew;
            
            A = handles.MarkersC.rt_down(Nbeats,ichan);
            Anew = MarkersC2.rt_down;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.rt_down(Nbeats,ichan) = Anew;
            
            A = handles.MarkersC.rt_up(Nbeats,ichan);
            Anew = MarkersC2.rt_up;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.rt_up(Nbeats,ichan) = Anew;
            
            
            w1 = max(w1,str2double(handles.tag_RTmin.String));
            w2 = min(w2,str2double(handles.tag_RTmax.String));
            
        case 6
            A = handles.MarkersC.tQRSoff(Nbeats,ichan);
            Anew = MarkersC2.tQRSoff;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.tQRSoff(Nbeats,ichan) = Anew;
            clear A Anew
            
            A = handles.MarkersC.tQRSon(Nbeats,ichan);
            Anew = MarkersC2.tQRSon;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.tQRSon(Nbeats,ichan) = Anew;
            clear A Anew
            
            A = handles.MarkersC.tQw(Nbeats,ichan);
            Anew = MarkersC2.tQw;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.tQw(Nbeats,ichan) = Anew;
            clear A Anew
            
            A = handles.MarkersC.tRw(Nbeats,ichan);
            Anew = MarkersC2.tRw;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.tRw(Nbeats,ichan) = Anew;
            T = Anew; % for plotting purposes
            clear A Anew
            
            A = handles.MarkersC.tSw(Nbeats,ichan);
            Anew = MarkersC2.tSw;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.tSw(Nbeats,ichan) = Anew;
            clear A Anew
            
            A = handles.MarkersC.QRS_amp(Nbeats,ichan);
            Anew = MarkersC2.QRS_amp;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.QRS_amp(Nbeats,ichan) = Anew;
            clear A Anew
            
            
            A = handles.MarkersC.QRS_area(Nbeats,ichan);
            Anew = MarkersC2.QRS_area;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.QRS_area(Nbeats,ichan) = Anew;
            clear A Anew
            
            A = handles.MarkersC.QRSoff(Nbeats,ichan);
            Anew = MarkersC2.QRSoff;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.QRSoff(Nbeats,ichan) = Anew;
            clear A Anew
            
            
            A = handles.MarkersC.QRSon(Nbeats,ichan);
            Anew = MarkersC2.QRSon;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.QRSon(Nbeats,ichan) = Anew;
            clear A Anew
            
            A = handles.MarkersC.QRSw(Nbeats,ichan);
            Anew = MarkersC2.QRSw;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.QRSw(Nbeats,ichan) = Anew;
            clear A Anew
            
            A = handles.MarkersC.QRSw_fw90(Nbeats,ichan);
            Anew = MarkersC2.QRSw_fw90;Anew(isnan(Anew)) = A(isnan(Anew)); % do not change if nan
            handles.MarkersC.QRSw_fw90(Nbeats,ichan) = Anew;
            clear A Anew
            
            w1 = max(w1,str2double(-handles.tag_DeltaT.String));
            w2 = min(w2,str2double(handles.tag_DTmax.String));
    end
end
clear A*
tag_Plot_Callback(hObject, [], handles);

% if handles.tag_all_beats.Value
% %     xs = find(xest>handles.spikes);xs = handles.spikes-handles.spikes(is);
%     for i = 1:length(Nbeats)
% %     plot(handles.tag_ax_sig,xs(i)+w1*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
% %     plot(handles.tag_ax_sig,xs(i)+w2*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
%     plot(handles.tag_ax_sig,handles.spikes(Nbeats(i))+w1*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
%     plot(handles.tag_ax_sig,handles.spikes(Nbeats(i))+w2*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
%
%     end
% else
%     plot(handles.tag_ax_sig,handles.spikes(Nbeats)+w1*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
%     plot(handles.tag_ax_sig,handles.spikes(Nbeats)+w2*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+mean(get(handles.tag_ax_sig,'ylim')),'--r')
% end
%

% normalize, just for visualization
if get(handles.tag_norm,'value')
    signal = signal./repmat(max(abs(signal)),[size(signal,1) 1]);
end
T = round(nanmedian(T,2)/1000*handles.ParamSig.frequency);
Am = zeros(size(Nbeats));
Am(~isnan(T)) = nanmedian(signal(T(~isnan(T)),ichan),2);

for i = 1:length(Nbeats)
    plot(handles.tag_ax_sig,handles.spikes(Nbeats(i))+w1*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+Am(i),'--r')
    plot(handles.tag_ax_sig,handles.spikes(Nbeats(i))+w2*[1 1],[-1 1]*range(get(handles.tag_ax_sig,'ylim'))/3+Am(i),'--r')
end

%%
% Update handles structure
guidata(hObject, handles);



% Autocorrection is based on:
% 1) Outliers are localized. They corresponds to a beat-to-beat variability higher than a threshold given by the user (Thrt)
% 2) For every CL repol. times are calculated inside a window centered around "x0Ref_tot" and "Trange" wide.
%   X0Ref and Trange can be given by the user or they can be calculated internally.
%   When they are given by the user they can be a scalar or a vector. If
%   they are a scalar, the same value is used for all CLs. If they are
%   vecotors, they can be of the same length as the CL.
%   -> When "X0Ref" si computed internally, it is a vector of length
%   Nbeats, whose values are those of a line fitted trhough the data using
%   a robust method.
% 3) A new repol time is localized in the window. Note that if there are
% several candidates, and the absolute maximum of the first derivative is similar to
% other local maxima, the repol time will be based according to the minimal temporak distance from X0Ref
%
% do_RTw_correction = 0;
% do_RTa_correction = 0;
% do_DT_correction = 0;
%
% if get(handles.tag_show_DT,'value')
%     do_DT_correction = 1;
% else
%     do_DT_correction = 0;
% end
%
% if get(handles.tag_show_RT,'value')|get(handles.tag_show_Tpeak,'value')
%     rtrypes = get(handles.tag_method,'string');
%     rtrypes = rtrypes{get(handles.tag_method,'value')};
%     if isequal(rtrypes,'Wyatt')|isequal(rtrypes,'Method')
%         do_RTw_correction = 1;
%     else
%         do_RTa_correction = 1;
%     end
% else
%     do_RTw_correction = 0;
%     do_RTa_correction = 0;
% end
%
%
% %%  AUTOCORRECTION
% do_median = 1;
%
% if (do_RTw_correction|do_RTa_correction)&do_DT_correction
%     answ = inputdlg({'CL [ms]','N beats for median filter','Reference Range [ms]','Threshold for outliers (beat-to-beat differences higher than) [ms]'},'Parameters for RT & AT cmodification',1,{'600,550,500,450,400,350','20','60','15'});
% elseif (do_RTw_correction|do_RTa_correction)
%     answ = inputdlg({'CL [ms]','N beats for median filter','Reference Range [ms]','Threshold for outliers (beat-to-beat differences higher than) [ms]'},'Parameters for RT modification',1,{'600,550,500,450,400,350','20','60','15'});
% elseif do_DT_correction
%     answdt = inputdlg({'CL [ms]','Reference Range [ms]','Threshold for outliers (beat-to-beat differences higher than) [ms]'},'Parameters for AT modification',1,{'600,550,500,450,400,350','30','8'});
%     answ{1} = answdt{1};answ{2} = [];answ(3:4) = answdt(2:3);
% end
% %%
% if isempty(answ)
%     return
% end
% %%
% Nmed_filt = str2double(answ{2});
% if rem(Nmed_filt,2)~=0
%     Nmed_filt=Nmed_filt+1;
% end
% Thrt = str2double(answ{4}); % ms
% DTmax = round(str2double(get(handles.tag_DTmax,'string'))/1000*handles.ParamSig.frequency);
% min_RT= round(str2double(get(handles.tag_RTmin,'string'))/1000*handles.ParamSig.frequency);
% max_RT = round(str2double(get(handles.tag_RTmax,'string'))/1000*handles.ParamSig.frequency);
% %% this is just to consider that in few cases activation can occur before the spike (due to filtering)
% % Note that this does not introduce any delay, because the temporal position of each marker does not change
% DeltaT = str2double(get(handles.tag_DeltaT,'string'))/1000*handles.ParamSig.frequency; % ms
% max_RT = max_RT + DeltaT;
% min_RT = min_RT + DeltaT;
% DTmax = DTmax + DeltaT;
% Hend = 3*DeltaT;
%
% %% spikes are modified here to consider in X an interval longer then CL. Use handles.spikes if true spike position is needed
% spikes = handles.spikes(:)-DeltaT; %samples
% spikes_samp = round( (handles.spikes(:)-DeltaT)/1000*handles.ParamSig.frequency); %samples
% signal =  handles.signals_proc;
% dt = handles.MarkersC.dt;
% dti = dt-repmat(spikes(:),[1 size(dt,2)]);
% rtiW = handles.MarkersC.rt_Wyatt-repmat(spikes(:),[1 size(dt,2)]);
% rtidown = handles.MarkersC.rt_down-repmat(spikes(:),[1 size(dt,2)]);
% tTp = handles.MarkersC.tTpeak-repmat(spikes(:),[1 size(dt,2)]);
% tTe = handles.MarkersC.tTend-repmat(spikes(:),[1 size(dt,2)]);
%
% %%
% hw = waitbar(0,'Initializing ...');
% CLtot = [];
% iicoma = find(answ{1}==',');
% if ~isempty(iicoma)&~isempty(answ{1})
%     CLtot(1) = str2double(answ{1}(1:iicoma(1)-1));
%     if length(iicoma)>1
%         for ic = 1 : length(iicoma)-1
%             CLtot(ic+1) = str2double(answ{1}(iicoma(ic)+1 : iicoma(ic+1)-1));
%         end
%         CLtot(ic+2) = str2double(answ{1}(iicoma(ic+1)+1 : end));
%     else
%         CLtot(2) = str2double(answ{1}(iicoma(1)+1 : end));
%     end
% else
%     CLtot = str2double(answ{1});
% end
% clear iicoma
%
% %% Define range around reference values
% Trange_tot = [];
% iicoma = find(answ{3}==',');
% if ~isempty(iicoma)
%     Trange_tot(1) = str2double(answ{3}(1:iicoma(1)-1));
%     if length(iicoma)>1
%         for ic = 1 : length(iicoma)-1
%             Trange_tot(ic+1) = str2double(answ{3}(iicoma(ic)+1 : iicoma(ic+1)-1));
%         end
%         Trange_tot(ic+2) = str2double(answ{3}(iicoma(ic+1)+1 : end));
%     else
%         Trange_tot(2) = str2double(answ{3}(iicoma(1)+1 : end));
%     end
% else
%     Trange_tot = str2double(answ{3});
% end
% clear iicoma
% %%
% if length(Trange_tot)>1&length(Trange_tot)~=length(CLtot)
%     error('Trange should be either a scalar or a vector with same dimensions as CL')
% end
%
% if length(CLtot)>1&length(Trange_tot)==1
%     Trange_tot = ones(1,length(CLtot))*Trange_tot;
% end
%
% %%
% for iCL = 1:length(CLtot)
%     CL = CLtot(iCL);
%     Trange = round(Trange_tot(iCL)/1000*handles.ParamSig.frequency); % samples
%
%     if ~isnan(CL)
%         Nbeats = find(diff(spikes)<CL+10 & diff(spikes)>CL-10);
%         if isempty(Nbeats)
%             hw2 = warndlg(['No heart beats with CL=',num2str(CL),'ms']);
%             pause(1)
%             close(hw2)
%             continue
%         else
%             Nbeats = [Nbeats;Nbeats(end)+1];
%         end
%     else
%         Nbeats = [1:length(spikes)];
%     end
%
%     iichan = get(handles.tag_List,'value');
%     Vdti = dti(Nbeats(2:end),:)-dti(Nbeats(1:end-1),:); % ms
%     Vrtiup = rtiW(Nbeats(2:end),:)-rtiW(Nbeats(1:end-1),:); % ms
%     Vrtidow = rtidown(Nbeats(2:end),:)-rtidown(Nbeats(1:end-1),:); % ms
%     Vtp = tTp(Nbeats(2:end),:)-tTp(Nbeats(1:end-1),:); % ms
%     Vte = tTe(Nbeats(2:end),:)-tTe(Nbeats(1:end-1),:); % ms
%     %% Create matrix [time,heart beat,electrode] (just to save time)
%     maxCL = round(800/1000*handles.ParamSig.frequency);
%     L=min([max(diff(spikes_samp)),maxCL]) + Hend;
%     X = nan(L-1,length(spikes_samp),length(iichan));
%     for i=1:length(spikes_samp)-1
%         H = spikes_samp(i):spikes_samp(i+1)-1+ Hend;
%         H(L:end)=[];H(H<1)=[];H(H>size(signal,1))=[];
%         X(1:length(H),i,:) = signal(H,iichan);
%     end
%     H = spikes_samp(i+1):size(signal,1);
%     H(L:end)=[];
%     X(1:length(H),i+1,:) = signal(H,iichan);
%     Xd = diff(X); Xd2 = diff(Xd); Xd3 = diff(Xd2);% Derivatives
%
%     %% REPOLARIZATION (up) for both Wyatt and Alternative
%     if do_RTw_correction|do_RTa_correction
%         clear ii j*
%         W = [min_RT:max_RT]; % windows in which rep time can be estimated
%         WTe = [min_RT : max_RT+round(100/1000*handles.ParamSig.frequency)];
%         W(W>(size(Xd,1)-2))=[];
%         WTe(WTe>(size(Xd2,1)-2))=[];
%
%         ii = Xd2(W(1:end-1),:,:).*Xd2(W(2:end),:,:)<0 & Xd(W(2:end),:,:)>0 & Xd3(W(2:end),:,:)<0; % zeros of II derivative & I derivative positive and II derivative negative (dV/dT max)
%         iiMax = Xd(WTe(1:end-1),:,:).*Xd(WTe(2:end),:,:)<0 & Xd2(WTe(2:end),:,:)<0; % zeros of I derivative & II derivative negative (V max)
%         iiMin = Xd(WTe(1:end-1),:,:).*Xd(WTe(2:end),:,:)<0 & Xd2(WTe(2:end),:,:)>0; % zeros of I derivative & II derivative positive (V min)
%
%         %% UP
%         for jchan = 1:length(iichan)
%             waitbar(jchan/length(iichan),hw,['[Correcting Repolarization (up) [CL=',num2str(CL),' #',num2str(iichan(jchan)),']...']);
%
%             rt = rtiW(Nbeats,iichan(jchan));
%             dt = dti(Nbeats,iichan(jchan));
%             if mean(isnan(rt))>0.3 | mean(isnan(dt))>0.3
%                 continue
%             end
%             % find outliers: values whose beat-to-beat change is higher than Thrt
%             % [samples]
%             iibeats = find(abs(Vrtiup(:,iichan(jchan)))>Thrt)+1;
%             iibeats = unique([iibeats(:)-1;iibeats(:);iibeats(:)+1;find(isnan(rt))]);
%             iibeats(iibeats>length(Nbeats)|iibeats==0|iibeats>length(rt))=[];
%             if do_median
%                 if length(rt)<=30
%                     xref = nanmedian(rt)*ones(1,length(rt));
%                 else
%                     x = [1:length(rt)];
%                     if sum(~isnan(rt))>0
%                         rt(1) = nanmedian(rt(1:5));rt(end) = nanmedian(rt(end-5:end));
%                         rt = interp1(x(~isnan(rt)),rt(~isnan(rt)),x);
%                     end
%                     xref = medfilt1([rt(1:round(Nmed_filt)/2) rt rt(end-round(Nmed_filt)/2:end)],Nmed_filt);
%                     xref = xref(round(Nmed_filt)/2+1:length(rt)+round(Nmed_filt)/2);
%                 end
%             else
%                 %             elseif sum(isnan(rt))<5&sum(length(rt))>5
%                 %                 [cc,gg] = fit([1:length(rt)]',rt(:),'poly1','exclude',isnan(rt),'robust','LAR');
%                 %                 xref =  cc.p1*[1:length(rt)] + cc.p2;
%                 %             elseif sum(isnan(rt))>5&sum(length(rt))>5
%                 %                 display(['Too many nans in : CL=',num2str(CL),' #',num2str(iichan(jchan))])
%                 %                 continue
%                 %             else
%                 %                 xref = nanmedian(rt)*ones(1,length(rt)); % ms
%             end
%
%
%             for j = 1:length(iibeats)
%                 tt = find(ii(:,Nbeats(iibeats(j)),jchan))+W(1)-1; % samples
%                 %                 tt(tt> xref(iibeats(j))/1000*handles.ParamSig.frequency + Trange/1000*handles.ParamSig.frequency/2 | tt<xref(iibeats(j))/1000*handles.ParamSig.frequency - Trange/2/1000*handles.ParamSig.frequency) = []; % Only consider a window centered around xest and WinMod wide
%                 tt(tt> xref(iibeats(j))/1000*handles.ParamSig.frequency + Trange/2 | tt<xref(iibeats(j))/1000*handles.ParamSig.frequency - Trange/2) = []; % Only consider a window centered around xest and WinMod wide
%                 %                 tt(tt>dt(iibeats(j))+max_ARI | tt<dt(iibeats(j))+min_ARI) = [];
%                 tt(tt>max_RT | tt<min_RT) = [];
%                 if ~isempty(tt)
%                     if length(tt)==1
%                         rtiW(Nbeats(iibeats(j)),iichan(jchan)) = tt/handles.ParamSig.frequency*1000;
%                     else
%                         Ratio = sort(Xd(tt,Nbeats(iibeats(j)),jchan)/max(Xd(tt,Nbeats(iibeats(j)),jchan)),'descend');
%                         if Ratio(2)<0.9
%                             [~,kk]=max(Xd(tt,Nbeats(iibeats(j)),jchan));
%                             tt = tt(kk);
%                             rtiW(Nbeats(iibeats(j)),iichan(jchan)) = tt/handles.ParamSig.frequency*1000;
%                         else
%                             [~,kk]= min(abs(tt-xref(iibeats(j))/1000*handles.ParamSig.frequency));
%                             tt = tt(kk);
%                             rtiW(Nbeats(iibeats(j)),iichan(jchan)) = tt/handles.ParamSig.frequency*1000;
%                         end
%                     end
%
%                     if ~handles.MarkersC.iiTwpos(iichan(jchan))
%                         % Update Tend
%                         x0 = tt;
%                         p1up=Xd(tt,Nbeats(iibeats(j)),jchan);
%                         tte = WTe(find(iiMax(:,Nbeats(iibeats(j)),jchan)))+1;
%                         tte(tte<x0)=[];
%                         if j==length(spikes_samp)
%                             tte2 = min(spikes_samp(j)+nanmedian(spikes_samp)-200/1000*handles.ParamSig.frequency -spikes_samp(j),size(Xd2,1));
%                         else
%                             tte2 = spikes_samp(j+1)-200/1000*handles.ParamSig.frequency -spikes_samp(j);
%                         end
%                         Tm = tte:tte2;
%                         Tm(abs(Xd(tte:tte2,Nbeats(iibeats(j)),jchan))>abs(0.1*p1up))=[];
%                         it = find(abs(diff(Tm))>1);
%                         if ~isempty(it)
%                             Tm(it(1)+1:end)=[];
%                         end
%                         V0 = nanmean(X(Tm,Nbeats(iibeats(j)),jchan));
%
%                         if isnan(V0); V0=0; end
%
%
%                         r = [-(X(tt,Nbeats(iibeats(j)),jchan)-p1up*tt-V0)/p1up];%r = roots([p1up (X(tt,Nbeats(iibeats(j)),jchan)-p1up*tt)]);
%                         if r<size(X,1) % no longer than cycle length
%                             tTe(Nbeats(iibeats(j)),iichan(jchan)) = r/handles.ParamSig.frequency*1000; % time for which the line tangent to tMin is equal to zero (for positive T waves)
%                         end
%                         clear r p1up
%                     end
%                 end
%             end
%             clear xref
%         end
%
%         %% DOWN (Alternative)
%         clear rt
%         for jchan = 1:length(iichan)
%             waitbar(jchan/length(iichan),hw,['[Correcting Repolarization (down) [CL=',num2str(CL),' #',num2str(iichan(jchan)),']...']);
%
%             rt = rtidown(Nbeats,iichan(jchan));
%             dt = dti(Nbeats,iichan(jchan));
%             if mean(isnan(rt))>0.3 | mean(isnan(dt))>0.3
%                 continue
%             end
%             % find outliers: values whose beat-to-beat change is higher than Thrt
%             % [samples]
%             iibeats = find(abs(Vrtidow(:,iichan(jchan)))>Thrt)+1;
%             iibeats = unique([iibeats(:)-1;iibeats(:);iibeats(:)+1;find(isnan(rt))]);
%             iibeats(iibeats>length(Nbeats)|iibeats==0|iibeats>length(rt))=[];
%             if do_median
%                 if length(rt)<=30
%                     xref = nanmedian(rt)*ones(1,length(rt));
%                 else
%                     x = [1:length(rt)];
%                     if sum(~isnan(rt))>0
%                         rt(1) = nanmedian(rt(1:5));rt(end) = nanmedian(rt(end-5:end));
%                         rt = interp1(x(~isnan(rt)),rt(~isnan(rt)),x);
%                     end
%                     xref = medfilt1([rt(1:round(Nmed_filt)/2) rt rt(end-round(Nmed_filt)/2:end)],Nmed_filt);
%                     xref = xref(round(Nmed_filt)/2+1:length(rt)+round(Nmed_filt)/2);
%                 end
%                 %             elseif sum(isnan(rt))<5&sum(length(rt))>5
%                 %                 [cc,gg] = fit([1:length(rt)]',rt(:),'poly1','exclude',isnan(rt),'robust','LAR');
%                 %                 xref =  cc.p1*[1:length(rt)] + cc.p2;
%                 %             elseif sum(isnan(rt))>5&sum(length(rt))>5
%                 %                 display(['Too many nans in : CL=',num2str(CL),' #',num2str(iichan(jchan))])
%                 %                 continue
%                 %             else
%                 %                 xref = nanmedian(rt)*ones(1,length(rt)); % ms
%             end
%
%             clear ii
%             ii_alt = Xd2(W(1:end-1),:,:).*Xd2(W(2:end),:,:)<0 & Xd(W(2:end),:,:)<0 & Xd3(W(2:end),:,:)>=0;; % zeros of II derivative & I derivative positive & III derivative > 0 (dV/dT min)
%
%             for j = 1:length(iibeats)
%                 tt = find(ii_alt(:,Nbeats(iibeats(j)),jchan))+W(1)-1; % samples
%                 %                 tt(tt> xref(iibeats(j))/1000*handles.ParamSig.frequency + Trange/1000*handles.ParamSig.frequency/2 | tt<xref(iibeats(j))/1000*handles.ParamSig.frequency - Trange/2/1000*handles.ParamSig.frequency) = []; % Only consider a window centered around xest and WinMod wide
%                 tt(tt> xref(iibeats(j))/1000*handles.ParamSig.frequency + Trange/2 | tt<xref(iibeats(j))/1000*handles.ParamSig.frequency - Trange/2) = []; % Only consider a window centered around xest and WinMod wide
%                 tt(tt>dt(iibeats(j))+max_ARI | tt<dt(iibeats(j))+min_ARI) = [];
%                 if ~isempty(tt)
%                     if length(tt)==1
%                         rtidown(Nbeats(iibeats(j)),iichan(jchan)) = tt/handles.ParamSig.frequency*1000;
%                     else
%                         Ratio = sort(Xd(tt,Nbeats(iibeats(j)),jchan)/max(Xd(tt,Nbeats(iibeats(j)),jchan)),'descend');
%                         if Ratio(2)<0.9
%                             [~,kk]=max(Xd(tt,Nbeats(iibeats(j)),jchan));
%                             tt = tt(kk);
%                             rtidown(Nbeats(iibeats(j)),iichan(jchan)) = tt/handles.ParamSig.frequency*1000;
%                         else
%                             [~,kk]= min(abs(tt-xref(iibeats(j))/1000*handles.ParamSig.frequency));
%                             tt = tt(kk);
%                             rtidown(Nbeats(iibeats(j)),iichan(jchan)) = tt/handles.ParamSig.frequency*1000;
%                         end
%                     end
%
%                     % Update Tend
%                     if handles.MarkersC.iiTwpos(iichan(jchan))
%                         x0 = tt;
%                         p1down=Xd(tt,Nbeats(iibeats(j)),jchan);
%                         tte = WTe(find(iiMin(:,Nbeats(iibeats(j)),jchan)))+1;
%                         tte(tte<x0)=[];
%                         if j==length(spikes_samp)
%                             tte2 = min(spikes_samp(j)+nanmedian(spikes_samp)-200/1000*handles.ParamSig.frequency -spikes_samp(j),size(Xd2,1));
%                         else
%                             tte2 = spikes_samp(j+1)-200/1000*handles.ParamSig.frequency -spikes_samp(j);
%                         end
%                         Tm = tte:tte2;
%                         Tm(abs(Xd(tte:tte2,Nbeats(iibeats(j)),jchan))>abs(0.1*p1down))=[];
%                         it = find(abs(diff(Tm))>1);
%                         if ~isempty(it)
%                             Tm(it(1)+1:end)=[];
%                         end
%                         V0 = nanmean(X(Tm,Nbeats(iibeats(j)),jchan));
%                         r = [-(X(tt,Nbeats(iibeats(j)),jchan)-p1down*tt-V0)/p1down];%roots([p1down (X(tt,Nbeats(iibeats(j)),jchan)-p1down*tt)]);
%                         if r<size(X,1) % no longer than cycle length
%                             tTe(Nbeats(iibeats(j)),iichan(jchan)) = r/handles.ParamSig.frequency*1000; % time for which the line tangent to tMin is equal to zero (for positive T waves)
%                         end
%                         clear p1down r
%                     end
%                 end
%             end
%             clear xref Ratio
%         end
%
%         %% Tpeak
%         clear rt
%         for jchan = 1:length(iichan)
%             waitbar(jchan/length(iichan),hw,['[Correcting T-peak [CL=',num2str(CL),' #',num2str(iichan(jchan)),']...']);
%
%             tTp_serie = tTp(Nbeats,iichan(jchan));
%             if mean(isnan(tTp_serie))>0.3 | mean(isnan(tTp_serie))>0.3
%                 continue
%             end
%             % find outliers: values whose beat-to-beat change is higher than Thrt
%             % [samples]
%             iibeats = find(abs(Vtp(:,iichan(jchan)))>Thrt)+1;
%             iibeats = unique([iibeats(:)-1;iibeats(:);iibeats(:)+1;find(isnan(tTp_serie))]);
%             iibeats(iibeats>length(Nbeats)|iibeats==0|iibeats>length(tTp_serie))=[];
%             if do_median
%                 if length(tTp_serie)<=30
%                     xref = nanmedian(tTp_serie)*ones(1,length(tTp_serie));
%                 else
%                     x = [1:length(tTp_serie)];
%                     if sum(~isnan(tTp_serie))>0
%                         tTp_serie(1) = nanmedian(tTp_serie(1:5));tTp_serie(end) = nanmedian(tTp_serie(end-5:end));
%                         tTp_serie = interp1(x(~isnan(tTp_serie)),tTp_serie(~isnan(tTp_serie)),x);
%                     end
%                     xref = medfilt1([tTp_serie(1:round(Nmed_filt)/2) tTp_serie tTp_serie(end-round(Nmed_filt)/2:end)],Nmed_filt);
%                     xref = xref(round(Nmed_filt)/2+1:length(tTp_serie)+round(Nmed_filt)/2);
%                 end
%             else
%                 %             elseif sum(isnan(tTp_serie))<5&sum(length(tTp_serie))>5
%                 %                 [cc,gg] = fit([1:length(tTp_serie)]',tTp_serie(:),'poly1','exclude',isnan(tTp_serie),'robust','LAR');
%                 %                 xref =  cc.p1*[1:length(tTp_serie)] + cc.p2;
%                 %             elseif sum(isnan(tTp_serie))>5&sum(length(tTp_serie))>5
%                 %                 display(['Too many nans in : CL=',num2str(CL),' #',num2str(iichan(jchan))])
%                 %                 continue
%                 %             else
%                 %                 xref = nanmedian(tTp_serie)*ones(1,length(tTp_serie)); % ms
%             end
%
%             clear ii
%             iiMax = Xd(W(1:end-1),:,:).*Xd(W(2:end),:,:)<0 & Xd2(W(2:end),:,:)<0; % zeros of I derivative & II derivative negative (V max)
%             iiMin = Xd(W(1:end-1),:,:).*Xd(W(2:end),:,:)<0 & Xd2(W(2:end),:,:)>0; % zeros of I derivative & II derivative positive (V min)
%             for j = 1:length(iibeats)
%                 if handles.MarkersC.iiTwpos(iichan(jchan))
%                     tt = W(find(iiMax(:,Nbeats(iibeats(j)),jchan)))+1; % sabples
%                     %                     tt(tt<=rtiW(Nbeats(iibeats(j)),iichan(jchan))/1000*handles.ParamSig.frequency |tt>=rtidown(Nbeats(iibeats(j)),iichan(jchan))/1000*handles.ParamSig.frequency )=[];
%                 else
%                     tt = W(find(iiMin(:,Nbeats(iibeats(j)),jchan)))+1;
%                     %                     tt(tt>= (rtidown(Nbeats(iibeats(j)),iichan(jchan))+rtiW(Nbeats(iibeats(j)),iichan(jchan)))/2/1000*handles.ParamSig.frequency )=[];
%                 end
%                 tt(tt> xref(iibeats(j))/1000*handles.ParamSig.frequency + Trange/2 | tt<xref(iibeats(j))/1000*handles.ParamSig.frequency - Trange/2) = []; % Only consider a window centered around xest and WinMod wide
%                 %                 figure,
%                 %                 plot(X(:,Nbeats(iibeats(j)),jchan)),hold on,plot(tt,X(tt,Nbeats(iibeats(j)),jchan),'or')
%                 %                 plot(xref(iibeats(j))/1000*handles.ParamSig.frequency + Trange/2 *[1 1],get(gca,'ylim'),'r')
%                 %                 plot(xref(iibeats(j))/1000*handles.ParamSig.frequency - Trange/2 *[1 1],get(gca,'ylim'),'r')
%
%                 if ~isempty(tt)
%                     if length(tt)==1
%                         tTp(Nbeats(iibeats(j)),iichan(jchan)) = tt/handles.ParamSig.frequency*1000;
%                     else
%                         [Ratio,is] = sort(abs(X(tt,Nbeats(iibeats(j)),jchan))/max(abs(X(tt,Nbeats(iibeats(j)),jchan))),'descend');
%                         if Ratio(2)<0.9 % if the highest peak is at least 10% highest than the second highest than take this one, otherwise take the closest to the reference
%                             tTp(Nbeats(iibeats(j)),iichan(jchan)) = tt(is(1))/handles.ParamSig.frequency*1000;
%                         else
%                             [~,kk]= min(abs(tt-xref(iibeats(j))/1000*handles.ParamSig.frequency));
%                             tTp(Nbeats(iibeats(j)),iichan(jchan)) = tt(kk)/handles.ParamSig.frequency*1000;
%
%                         end
%                     end
%                 end
%                 %%
%
%             end
%         end
%
%     end
%
%
%
%     %% depolarization
%     if do_DT_correction
%         clear ii j*
%         ii = Xd2([1:DTmax+80],:,:).*Xd2([2:DTmax+81],:,:)<0 & Xd([1:DTmax+80],:,:)<0 & Xd3([1:DTmax+80],:,:)>0; % zeros of II derivative & I derivative negative & III derivative > 0
%         for jchan = 1:length(iichan)
%
%             waitbar(jchan/length(iichan),hw,['[Correcting Activation (up) [CL=',num2str(CL),' #',num2str(iichan(jchan)),']...']);
%             dt = dti(Nbeats,iichan(jchan));
%
%             % find outliers: values whose beat-to-beat change is higher than Thrt
%             % [samples]
%             iibeats = find(abs(Vdti(:,iichan(jchan)))>Thrt)+1;
%             iibeats = unique([iibeats(:)-1;iibeats(:);iibeats(:)+1;find(isnan(dt))]);
%             iibeats(iibeats>length(Nbeats)|iibeats==0)=[];
%             xref = nanmedian(dt);
%
%             if sum(ii(:))>0
%                 %% Depolarization time
%                 for j = 1:length(iibeats) % j=heart beat
%                     tt = find(ii(:,Nbeats(iibeats(j)),jchan)); % samples
%                     tt(tt>xref/1000*handles.ParamSig.frequency+Trange/2|tt<xref/1000*handles.ParamSig.frequency-Trange/2) = []; % Only consider a window centered around xest and WinMod wide
%
%                     if ~isempty(tt)
%                         [~,kk] = min(Xd(tt,j,jchan)); % min inside the window defined by WinMod
%                         dti(Nbeats(iibeats(j)),iichan(jchan)) = tt(kk)/handles.ParamSig.frequency*1000;
%                     end
%                 end
%             end
%
%         end
%     end
%
%
% end
%
%
% handles.MarkersC.rt_up = rtiW+repmat(spikes,[1 size(rtiW,2)]);
% handles.MarkersC.rt_down = rtidown+repmat(spikes,[1 size(rtidown,2)]);
% handles.MarkersC.dt = dti+repmat(spikes,[1 size(rtidown,2)]);;
% handles.MarkersC.tTpeak = tTp+repmat(spikes,[1 size(tTp,2)]);
% handles.MarkersC.tTend = tTe+repmat(spikes,[1 size(tTp,2)]);
%
%
% iiTwpos = logical(handles.MarkersC.iiTwpos);
% handles.MarkersC.rt_Wyatt = handles.MarkersC.rt_up;
% handles.MarkersC.rt_Alternative(:,~iiTwpos) = handles.MarkersC.rt_up(:,~iiTwpos);
% handles.MarkersC.rt_Alternative(:,iiTwpos) = handles.MarkersC.rt_down(:,iiTwpos);
% close(hw)
%
% hold(handles.tag_ax_sig,'off')
% axes(handles.tag_ax_sig)
% plot(handles.tag_ax_sig,[0.4 0.6],[0 0])
% text(.5,.05,'Re-Plot')
% hold(handles.tag_ax_spikes,'off')
% axes(handles.tag_ax_spikes)
% plot(handles.tag_ax_spikes,[0.4 0.6],[0 0])
% text(.5,.05,'Re-Plot')
%
% %
% guidata(hObject,handles);
%
% tag_Plot_Callback(hObject,[], handles);

% --- Executes on button press in tag_show_DT.
function tag_show_DT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_DT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_DT


% --- Executes on button press in tag_show_ARI.
function tag_show_ARI_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_ARI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_ARI


% --- Executes on button press in tag_show_RT.
function tag_show_RT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_RT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_RT


% --- Executes on button press in tag_legend.
function tag_legend_Callback(hObject, eventdata, handles)
% hObject    handle to tag_legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_legend
cc0 = findobj(0,'name','Markers_GUI_v2');
if ~isempty(cc0)
    cc1 = findobj(cc0(1),'tag','legend');
    if handles.tag_legend.Value
        set(cc1,'visible','on')
    else
        set(cc1,'visible','off')
    end
end


% --- Executes on button press in tag_all_chan.
function tag_all_chan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_all_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_all_chan


% --- Executes on button press in tag_export_plot.
function tag_export_plot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_export_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



Fig2 = figure;
ax2(1)=copyobj(handles.tag_ax_sig, Fig2);
ax2(2)=copyobj(handles.tag_ax_spikes, Fig2);
set(ax2(2),'position',[.1 .1 .85 .38])
set(ax2(1),'position',[.1 .53 .85 .38])
%%
legend( get(ax2(2),'children'),'DT','ARI','RT','ISO','Tp','Tend')
iv = get(handles.tag_List,'value');
for j = 1:length(iv)
    ichan = iv(j);
    Leg(j) =  {[handles.ParamSig.Label{ichan},'(',num2str(handles.SNR(ichan),2),'dB)']};
end

legend(ax2(1),Leg,'AutoUpdate','off')


% --- Executes on button press in tag_norm.
function tag_norm_Callback(hObject, eventdata, handles)
% hObject    handle to tag_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_norm


% --- Executes on button press in tag_delete.
function tag_delete_Callback(hObject, eventdata, handles)
% hObject    handle to tag_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dcm_obj = datacursormode(handles.tag_figure_markers);
set(dcm_obj,'updateFcn',[],'enable','off')
info_struct = getCursorInfo(dcm_obj);
delete(findall(handles.tag_figure_markers,'Type','hggroup','HandleVisibility','off'));
doplot = 0;

if isempty(info_struct)
    return
end

for iter = 1:length(info_struct)
    x = info_struct(iter).Position(1);
    ix = round(info_struct(iter).Position(1)*handles.ParamSig.frequency/1000);
    
    
    %%
    % Filtered/not filtered
    if get(handles.tag_show_filter_RAW,'value')
        signal=handles.S;
    elseif get(handles.tag_show_filter_RT,'value')
        signal=handles.signals_proc;
    elseif get(handles.tag_show_filter_AT,'value')
        signal=handles.signals_proc_AT;
    else
        h = helpdlg('Please select one of the filter mode: RAW, RT, AT');
        return
    end
    
    
    % dV/dt
    if get(handles.tag_dV,'value')
        signal=diff(signal);
    end
    if get(handles.tag_norm,'value')
        signal = signal./repmat(max(abs(signal)),[size(signal,1) 1]);
    end
    %%
    iil = handles.tag_List.Value;
    ichan = find(signal(ix,iil)==info_struct(iter).Position(2));
    ichan = iil(ichan);
    
    do_all_chan = get(handles.tag_all_chan,'value');
    do_all_beats = get(handles.tag_all_beats,'value');
    do_all_beats_sameCL = get(handles.tag_ALL_beat_CL,'value');
    if do_all_beats
        answ = questdlg('Delete markers from all beats?');
        if ~isequal(answ,'Yes');
            return
        end
    end
    
    if do_all_beats_sameCL
        answ = questdlg('Delete markers from all beats [with same CL]?');
        if ~isequal(answ,'Yes');
            return
        end
    end
    
    ichan_all = get(handles.tag_List,'value');
    % - isot
    ii=find(handles.MarkersC.isot(:,ichan)==x);
    if do_all_beats&~isempty(ii)
        ii = 1:size(handles.MarkersC.isot,1);
        if do_all_chan
            ix = round(handles.MarkersC.isot(:,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.isot(:,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if do_all_beats_sameCL&~isempty(ii)
        cl = handles.spikes(ii)-handles.spikes(ii-1);
        ii = find(diff(handles.spikes)<cl+10 & diff(handles.spikes)>cl-10);
        ii = [ii(:);ii(end)+1];
        if do_all_chan
            ix = round(handles.MarkersC.isot(ii,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.isot(ii,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if ~isempty(ii)
        doplot = 1;
        if do_all_chan
            answ = questdlg('Delete iso-t from all visualised channels?');
            if isequal(answ,'Yes');
                xALL = handles.MarkersC.isot(ii,ichan_all);
                handles.MarkersC.isot(ii,ichan_all)=nan;
            else
                return
            end
        else
            handles.MarkersC.isot(ii,ichan)=nan;
        end
    end
    % - rtA
    ii=find(handles.MarkersC.rt_Alternative(:,ichan)==x);
    if do_all_beats&~isempty(ii)
        ii = 1:size(handles.MarkersC.rt_Alternative,1);
        if do_all_chan
            ix = round(handles.MarkersC.rt_Alternative(:,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.rt_Alternative(:,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if do_all_beats_sameCL&~isempty(ii)
        cl = handles.spikes(ii)-handles.spikes(ii-1);
        ii = find(diff(handles.spikes)<cl+10 & diff(handles.spikes)>cl-10);
        ii = [ii(:);ii(end)+1];
        if do_all_chan
            ix = round(handles.MarkersC.rt_Alternative(ii,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.rt_Alternative(ii,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if ~isempty(ii)
        doplot = 1;
        if do_all_chan
            answ = questdlg('Delete RT-A from all visualised channels?');
            if isequal(answ,'Yes');
                xALL = handles.MarkersC.rt_Alternative(ii,ichan_all);
                handles.MarkersC.rt_Alternative(ii,ichan_all)=nan;
            else
                return
            end
        else
            handles.MarkersC.rt_Alternative(ii,ichan)=nan;
        end
    end
    % - rtW
    ii=find(handles.MarkersC.rt_Wyatt(:,ichan)==x);
    if do_all_beats&~isempty(ii)
        ii = 1:size(handles.MarkersC.rt_Wyatt,1);
        if do_all_chan
            ix = round(handles.MarkersC.rt_Wyatt(:,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.rt_Wyatt(:,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if do_all_beats_sameCL&~isempty(ii)
        cl = handles.spikes(ii)-handles.spikes(ii-1);
        ii = find(diff(handles.spikes)<cl+10 & diff(handles.spikes)>cl-10);
        ii = [ii(:);ii(end)+1];
        if do_all_chan
            ix = round(handles.MarkersC.rt_Wyatt(ii,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.rt_Wyatt(ii,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if ~isempty(ii)
        doplot = 1;
        if do_all_chan
            answ = questdlg('Delete RT-W from all visualised channels?');
            if isequal(answ,'Yes');
                xALL = handles.MarkersC.rt_Wyatt(ii,ichan_all);
                handles.MarkersC.rt_Wyatt(ii,ichan_all)=nan;
            else
                return
            end
        else
            handles.MarkersC.rt_Wyatt(ii,ichan)=nan;
        end
    end
    % - rt_down
    ii=find(handles.MarkersC.rt_down(:,ichan)==x);
    if do_all_beats&~isempty(ii)
        ii = 1:size(handles.MarkersC.rt_down,1);
        if do_all_chan
            ix = round(handles.MarkersC.rt_down(:,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.rt_down(:,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if do_all_beats_sameCL&~isempty(ii)
        cl = handles.spikes(ii)-handles.spikes(ii-1);
        ii = find(diff(handles.spikes)<cl+10 & diff(handles.spikes)>cl-10);
        ii = [ii(:);ii(end)+1];
        if do_all_chan
            ix = round(handles.MarkersC.rt_down(ii,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.rt_down(ii,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if ~isempty(ii)
        doplot = 1;
        if do_all_chan
            answ = questdlg('Delete RT-down from all visualised channels?');
            if isequal(answ,'Yes');
                xALL = handles.MarkersC.rt_down(ii,ichan_all);
                handles.MarkersC.rt_down(ii,ichan_all)=nan;
            else
                return
            end
        else
            handles.MarkersC.rt_down(ii,ichan)=nan;
        end
    end
    % - rt_up
    ii=find(handles.MarkersC.rt_up(:,ichan)==x);
    if do_all_beats&~isempty(ii)
        ii = 1:size(handles.MarkersC.rt_up,1);
        if do_all_chan
            ix = round(handles.MarkersC.rt_up(:,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.rt_up(:,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if do_all_beats_sameCL&~isempty(ii)
        cl = handles.spikes(ii)-handles.spikes(ii-1);
        ii = find(diff(handles.spikes)<cl+10 & diff(handles.spikes)>cl-10);
        ii = [ii(:);ii(end)+1];
        if do_all_chan
            ix = round(handles.MarkersC.rt_up(ii,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.rt_up(ii,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if ~isempty(ii)
        doplot = 1;
        if do_all_chan
            answ = questdlg('Delete RT-down from all visualised channels?');
            if isequal(answ,'Yes');
                xALL = handles.MarkersC.rt_up(ii,ichan_all);
                handles.MarkersC.rt_up(ii,ichan_all)=nan;
            else
                return
            end
        else
            handles.MarkersC.rt_up(ii,ichan)=nan;
        end
    end
    % - dt
    ii=find(handles.MarkersC.dt(:,ichan)==x);
    if do_all_beats&~isempty(ii)
        ii = 1:size(handles.MarkersC.dt,1);
        if do_all_chan
            ix = round(handles.MarkersC.dt(:,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.dt(:,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if do_all_beats_sameCL&~isempty(ii)
        cl = handles.spikes(ii)-handles.spikes(ii-1);
        ii = find(diff(handles.spikes)<cl+10 & diff(handles.spikes)>cl-10);
        ii = [ii(:);ii(end)+1];
        if do_all_chan
            ix = round(handles.MarkersC.dt(ii,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.dt(ii,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if ~isempty(ii)
        doplot = 1;
        if do_all_chan
            answ = questdlg('Delete DT from all visualised channels?');
            if isequal(answ,'Yes');
                xALL = handles.MarkersC.dt(ii,ichan_all);
                handles.MarkersC.dt(ii,ichan_all)=nan;
            else
                return
            end
        else
            handles.MarkersC.dt(ii,ichan)=nan;
        end
    end
    % - tTpeak
    ii=find(handles.MarkersC.tTpeak(:,ichan)==x);
    if do_all_beats&~isempty(ii)
        ii = 1:size(handles.MarkersC.tTpeak,1);
        if do_all_chan
            ix = round(handles.MarkersC.tTpeak(:,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.tTpeak(:,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if do_all_beats_sameCL&~isempty(ii)
        cl = handles.spikes(ii)-handles.spikes(ii-1);
        ii = find(diff(handles.spikes)<cl+10 & diff(handles.spikes)>cl-10);
        ii = [ii(:);ii(end)+1];
        if do_all_chan
            ix = round(handles.MarkersC.tTpeak(ii,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.tTpeak(ii,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if ~isempty(ii)
        doplot = 1;
        if do_all_chan
            answ = questdlg('Delete Tpeak from all visualised channels?');
            if isequal(answ,'Yes');
                xALL = handles.MarkersC.tTpeak(ii,ichan_all);
                handles.MarkersC.tTpeak(ii,ichan_all)=nan;
            else
                return
            end
        else
            handles.MarkersC.tTpeak(ii,ichan)=nan;
        end
    end
    %%
    % - tTend
    ii=find(handles.MarkersC.tTend(:,ichan)==x);
    if do_all_beats&~isempty(ii)
        ii = 1:size(handles.MarkersC.tTend,1);
        if do_all_chan
            ix = round(handles.MarkersC.tTend(:,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.tTend(:,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if do_all_beats_sameCL&~isempty(ii)
        cl = handles.spikes(ii)-handles.spikes(ii-1);
        ii = find(diff(handles.spikes)<cl+10 & diff(handles.spikes)>cl-10);
        ii = [ii(:);ii(end)+1];
        if do_all_chan
            ix = round(handles.MarkersC.tTend(ii,ichan_all)/1000*handles.ParamSig.frequency);
        else
            ix = round(handles.MarkersC.tTend(ii,ichan)/1000*handles.ParamSig.frequency);
        end
    end
    if ~isempty(ii)
        doplot = 1;
        if do_all_chan
            answ = questdlg('Delete Tpeak from all visualised channels?');
            if isequal(answ,'Yes');
                xALL = handles.MarkersC.tTend(ii,ichan_all);
                handles.MarkersC.tTend(ii,ichan_all)=nan;
            else
                return
            end
        else
            handles.MarkersC.tTend(ii,ichan)=nan;
        end
    end
    %%
    
    if doplot
        hold(handles.tag_ax_sig,'on')
        if do_all_chan&(~do_all_beats)&(~do_all_beats_sameCL)
            ixALL = round(xALL*handles.ParamSig.frequency/1000);
            iinan = find(isnan(xALL));xALL(iinan) = []; ixALL(iinan) = [];
            
            a = nanmean(xALL);
            b = nanmean(nanmean(signal(ixALL,ichan_all)));
            plot(handles.tag_ax_sig,a,b,'r*','linewidth',4,'markersize',30)
        elseif do_all_chan&(do_all_beats|do_all_beats_sameCL)
            for jc = 1:size(ix,2)
                a = ix(:,jc);
                a(isnan(a)) = [];
                a(a>size(signal,1))=[];
                b = signal(a,ichan_all(:,jc));
                plot(a/handles.ParamSig.frequency*1000,b,'rx','linewidth',3,'markersize',15)
            end
        else %~do_all_chan&do_all_beats
            a = ix;
            a(isnan(a)) = [];
            b = signal(a,ichan);
            plot(a/handles.ParamSig.frequency*1000,b,'rx','linewidth',3,'markersize',15)
            
        end
    else
        h = helpdlg('No point found');
        pause(1)
        close(h)
        return
    end
end
%% update
guidata(hObject,handles);


% --- Executes on slider movement.
function tag_slider_Callback(hObject, eventdata, handles)
% hObject    handle to tag_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
rr = size(handles.signals_proc,1);
xx = get(handles.tag_ax_sig,'xlim');
% sstep = get(handles.tag_slider,'sliderStep');
% if sstep(1)*rr>(xx(2)-xx(1))
%     set(handles.tag_slider,'sliderStep',[(xx(2)-xx(1))/rr/4 (xx(2)-xx(1))/rr]);
% else
%     set(handles.tag_slider,'sliderStep',[.01 .1]);
% end
set(handles.tag_ax_sig,'xlim',get(hObject,'Value')*rr+[0 xx(2)-xx(1)]);
% set(handles.tag_slider,'sliderStep',[.01 .1]);

% --- Executes during object creation, after setting all properties.
function tag_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in tag_show_Tpeak.
function tag_show_Tpeak_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_Tpeak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_Tpeak


% --- Executes on button press in tag_show_Tiso.
function tag_show_Tiso_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_Tiso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_Tiso


% --- Executes on button press in tag_summary.
function tag_summary_Callback(hObject, eventdata, handles)
% hObject    handle to tag_summary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = get(handles.tag_summary_type,'string');
snrlim = str2double(get(handles.tag_SNRth_L,'string'));

if isequal(str{get(handles.tag_summary_type,'value')},'Quality')
    figure(10)
    ax(1)=subplot(311);
    plot(handles.sig_corr.*handles.SNR,'squarek--','markerfacecolor','k');
    ax(2)=subplot(312);
    plot(handles.SNR,'o--','markerfacecolor','b');
    ax(3)=subplot(313);
    plot(handles.sig_corr,'o--','markerfacecolor','b');
    linkaxes(ax,'x')
    xlabel(ax(3),'Channels')
    ylabel(ax(1),['[dB*n.u.]']),ylabel(ax(2),['[dB]']),ylabel(ax(3),['[n.u.]'])
    title(ax(1),['SNR * Morphological Stability']),title(ax(2),['SNR']),title(ax(3),['Morphological Stability'])
elseif isequal(str{get(handles.tag_summary_type,'value')},'Markers(STV)')
    
    spdiff = abs(diff(diff(handles.spikes)));
    ibko = find(spdiff>0.010*handles.ParamSig.frequency);
    ibko = unique([ibko(:)+1;ibko(:)+2;ibko(:)+3]);
    ibko(ibko<1|ibko>length(handles.spikes))=[];
    
    dti=handles.MarkersC.dt-repmat(handles.spikes(:),[1,size(handles.MarkersC.dt,2)]);
    DTvar = nanmean(abs(diff(dti)))/handles.ParamSig.frequency*1000;
    
    rti=handles.MarkersC.rt_Wyatt-repmat(handles.spikes(:),[1,size(handles.MarkersC.dt,2)]);
    RTvar = nanmean(abs(diff(rti)))/handles.ParamSig.frequency*1000;
    
    figure(10)
    ax(1)=subplot(211);
    hold off
    aa(1,1)=plot(DTvar,'o--','markerfacecolor','b');hold on
    ii = find(DTvar(:)>5&handles.SNR(:)>snrlim&handles.sig_corr(:,1)>0.85);
    if ~isempty(ii)
        aa(2,1)=plot(ii,DTvar(ii),'o','markerfacecolor','r');
    else
        aa(2,1)=plot(nan,nan,'o','markerfacecolor','r');
    end
    kk = find(handles.SNR(:)<=snrlim|handles.sig_corr(:,1)<=0.85);
    if ~isempty(kk)
        aa(3,1)=plot(kk,DTvar(kk),'o','markerfacecolor',[.6 .6 .6],'markeredgecolor','k');
    end
    plot([0 (length(DTvar)+1)],[1 1]*5,'--k')
    ax(2)=subplot(212);
    hold off
    aa(1,2)=plot(RTvar,'o--','markerfacecolor','b');hold on
    ii = find(RTvar(:)>10&handles.SNR(:)>snrlim&handles.sig_corr(:,1)>0.85);
    if ~isempty(ii)
        aa(2,2)=plot(ii,RTvar(ii),'o','markerfacecolor','r');
    else
        aa(2,2)=plot(nan,nan,'o','markerfacecolor','r');
    end
    if ~isempty(kk)
        aa(3,2)=plot(kk,RTvar(kk),'o','markerfacecolor',[.6 .6 .6],'markeredgecolor','k');
    end
    plot([0 (length(RTvar)+1)],[1 1]*10,'--k')
    
    if ~isempty(kk)
        ll1 = legend([aa(2:3,1)],'DTV>5ms',['SNR<',num2str(snrlim),' dB  |corr<0.85']);
        ll2 = legend([aa(2:3,2)],'RTV>10ms',['SNR<',num2str(snrlim),' dB  |corr<0.85']);
    else
        ll1 = legend([aa(:,1)],'DTV>5ms',['SNR<',num2str(snrlim),' dB  |corr<0.85']);
        ll2 = legend([aa(:,2)],'RTV>10ms',['SNR<',num2str(snrlim),' dB  |corr<0.85']);
    end
    
    linkaxes(ax,'x')
    xlabel(ax(2),'Channels')
    ylabel(ax(1),['[ms]']),ylabel(ax(2),['[ms]'])
    title(ax(1),'Short term variability (ACTIVATION)'),title(ax(2),'Short term variability (REPOLARISATION)')
    
elseif isequal(str{get(handles.tag_summary_type,'value')},'DT(med)')
    dt = handles.MarkersC.dt - repmat(handles.spikes(:),[1,size(handles.MarkersC.dt,2)]);
    dtmed = nanmedian(dt,1);
    figure
    plot(dtmed,'o--');
    hold on
    aa=plot(find(handles.SNR<snrlim),dtmed(find(handles.SNR<snrlim)),'marker','o','markerfacecolor',[.8 .8 .8],'markeredgecolor',[.8 .8 .8],'linestyle','none');
    plot(find(handles.SNR>=snrlim),dtmed(find(handles.SNR>=snrlim)),'marker','o','markerfacecolor',[0 0 1],'linestyle','none')
    %     plot(dtmed,'o-'),
    xlabel('Channels'),ylabel(['[ms]']),title('Median Activartion Time')
    legend(aa,['SNR<',num2str(snrlim),'dB'])
elseif isequal(str{get(handles.tag_summary_type,'value')},'RT(med)')
    rtW = handles.MarkersC.rt_Wyatt - repmat(handles.spikes(:),[1,size(handles.MarkersC.dt,2)]);
    rtWmed = nanmedian(rtW,1);
    
    rtA = handles.MarkersC.rt_Alternative - repmat(handles.spikes(:),[1,size(handles.MarkersC.dt,2)]);
    rtAmed = nanmedian(rtA,1);
    
    figure
    ax(1) = subplot(211);
    %     plot(rtWmed,'o-'),
    plot(rtWmed,'o--');
    hold on
    aa=plot(find(handles.SNR<snrlim),rtWmed(find(handles.SNR<snrlim)),'marker','o','markerfacecolor',[.8 .8 .8],'markeredgecolor',[.8 .8 .8],'linestyle','none');
    plot(find(handles.SNR>=snrlim),rtWmed(find(handles.SNR>=snrlim)),'marker','o','markerfacecolor',[0 0 1],'linestyle','none')
    xlabel('Channels'),ylabel(['[ms]']),title('Median Repolarization (WYATT) Time')
    
    ax(2) = subplot(212);
    %     plot(rtAmed,'o-'),
    plot(rtAmed,'o--');
    hold on
    aa=plot(find(handles.SNR<snrlim),rtAmed(find(handles.SNR<snrlim)),'marker','o','markerfacecolor',[.8 .8 .8],'markeredgecolor',[.8 .8 .8],'linestyle','none');
    plot(find(handles.SNR>=snrlim),rtAmed(find(handles.SNR>=snrlim)),'marker','o','markerfacecolor',[0 0 1],'linestyle','none')
    xlabel('Channels'),ylabel(['[ms]']),title('Median Repolarization (ALTERNATIVE) Time')
    legend(aa,['SNR<',num2str(snrlim),'dB'])
    
elseif isequal(str{get(handles.tag_summary_type,'value')},'ARI(med)')
    ariW = handles.MarkersC.rt_Wyatt - handles.MarkersC.dt;
    ariWmed = nanmedian(ariW,1);
    
    ariA = handles.MarkersC.rt_Alternative - handles.MarkersC.dt;
    ariAmed = nanmedian(ariA,1);
    
    figure
    ax(1) = subplot(211);
    %     plot(ariWmed,'o-'),
    plot(ariWmed,'o--');
    hold on
    aa=plot(find(handles.SNR<snrlim),ariWmed(find(handles.SNR<snrlim)),'marker','o','markerfacecolor',[.8 .8 .8],'markeredgecolor',[.8 .8 .8],'linestyle','none');
    plot(find(handles.SNR>=snrlim),ariWmed(find(handles.SNR>=snrlim)),'marker','o','markerfacecolor',[0 0 1],'linestyle','none')
    xlabel('Channels'),ylabel(['[ms]']),title('Median ARI (WYATT) Time')
    
    ax(2) = subplot(212);
    %     plot(ariAmed,'o-'),
    plot(ariAmed,'o--');
    hold on
    aa=plot(find(handles.SNR<snrlim),ariAmed(find(handles.SNR<snrlim)),'marker','o','markerfacecolor',[.8 .8 .8],'markeredgecolor',[.8 .8 .8],'linestyle','none');
    plot(find(handles.SNR>=snrlim),ariAmed(find(handles.SNR>=snrlim)),'marker','o','markerfacecolor',[0 0 1],'linestyle','none')
    xlabel('Channels'),ylabel(['[ms]']),title('Median ARI (ALTERNATIVE) Time')
    legend(aa,['SNR<',num2str(snrlim),'dB'])
    
end









% --- Executes on selection change in tag_summary_type.
function tag_summary_type_Callback(hObject, eventdata, handles)
% hObject    handle to tag_summary_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_summary_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_summary_type


% --- Executes during object creation, after setting all properties.
function tag_summary_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_summary_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_ALL_beat_CL.
function tag_ALL_beat_CL_Callback(hObject, eventdata, handles)
% hObject    handle to tag_ALL_beat_CL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_ALL_beat_CL



function tag_W_mod_Callback(hObject, eventdata, handles)
% hObject    handle to tag_W_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_W_mod as text
%        str2double(get(hObject,'String')) returns contents of tag_W_mod as a double


% --- Executes during object creation, after setting all properties.
function tag_W_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_W_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_show_Tend.
function tag_show_Tend_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_Tend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_Tend



function tag_SNRth_U_Callback(hObject, eventdata, handles)
% hObject    handle to tag_SNRth_U (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_SNRth_U as text
%        str2double(get(hObject,'String')) returns contents of tag_SNRth_U as a double


% --- Executes during object creation, after setting all properties.
function tag_SNRth_U_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_SNRth_U (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end











% --- Executes on button press in tag_info_point.
function tag_info_point_Callback(hObject, eventdata, handles)
% hObject    handle to tag_info_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


dcm_obj = datacursormode(handles.tag_figure_markers);
info_struct = getCursorInfo(dcm_obj);

pos = info_struct.Position;
% ax = get(get(event_obj,'target'),'parent');
ax = get(info_struct.Target,'parent');

% handles = guidata(gcf);
iic = get(handles.tag_List,'value');


if ax == handles.tag_ax_sig
    
    xx = pos(1)-handles.spikes;xx(xx<0)=[];
    [~,im] = min(xx);
    tt = (pos(1)-handles.spikes(im));
    ip = round(pos(1)/1000*handles.ParamSig.frequency);
    if get(handles.tag_show_filter_RT,'value')
        if get(handles.tag_dV,'value')
            signal=diff(handles.signals_proc);
        else
            signal=handles.signals_proc;
        end
        if get(handles.tag_norm,'value')
            signal = signal(ip,iic)./max(abs(signal(:,iic)));
            jj = find(signal==pos(2));
        else
            jj = find(signal(ip,iic)==pos(2));
        end
    elseif get(handles.tag_show_filter_RAW,'value')
        if get(handles.tag_dV,'value')
            signal=diff(handles.S);
        else
            signal=handles.S;
        end
        if get(handles.tag_norm,'value')
            signal = signal(ip,iic)./max(abs(signal(:,iic)));
            jj = find(signal==pos(2));
        else
            jj = find(signal(ip,iic)==pos(2));
        end
    elseif get(handles.tag_show_filter_AT,'value')
        if get(handles.tag_dV,'value')
            signal=diff(handles.signals_proc_AT);
        else
            signal=handles.signals_proc_AT;
        end
        if get(handles.tag_norm,'value')
            signal = signal(ip,iic)./max(abs(signal(:,iic)));
            jj = find(signal==pos(2));
        else
            jj = find(signal(ip,iic)==pos(2));
        end
    end
    
    polarity = 'N/A';
    if isfield(handles,'MarkersC')
        if handles.MarkersC.iiTwpos(iic(jj(1)))==1
            polarity = 'pos';
        else
            polarity = 'neg';
        end
    else
        if isfield(handles,'Markers')
            if handles.Markers.iiTwpos(iic(jj(1)))==1
                polarity = 'pos';
            else
                polarity = 'neg';
            end
        end
    end
    
    if ~isempty(jj)
        output_txt= { ['Lab:',handles.ParamSig.Label{iic(jj(1))}],
            ['IC',num2str(iic(jj(1)))]
            ['t =',num2str(tt,3),' ms']
            ['Beat# = ',num2str(im)]
            ['x =',num2str(pos(1)),' samp']
            [polarity,' TW']};
        set(handles.tag_info_point_txt,'string',output_txt)
    end
    clear polarity
    
    
elseif ax == handles.tag_ax_spikes
    
    if get(handles.tag_show_RT,'value')&(get(handles.tag_method,'value')==2)
        rti = handles.MarkersC.rt_Wyatt(:,iic) - repmat(handles.spikes(:),[1 length(iic)]);
        [ii,jj] = find(handles.MarkersC.rt_Wyatt(:,iic)==pos(1)&rti == pos(2));
        if ~isempty(ii)
            output_txt= { ['Lab:',handles.ParamSig.Label{iic(jj(1))}],
                ['IC',num2str(iic(jj(1)))],
                ['RT=',num2str(rti(ii(1),jj(1)),3),'ms']
                };
            set(handles.tag_info_point_txt,'string',output_txt)
        end
    end
    
    if get(handles.tag_show_RT,'value')&(get(handles.tag_method,'value')==3)
        rti = handles.MarkersC.rt_Alternative(:,iic) - repmat(handles.spikes(:),[1 length(iic)]);
        [ii,jj] = find(handles.MarkersC.rt_Alternative(:,iic)==pos(1)&rti == pos(2));
        
        if ~isempty(ii)
            output_txt= { ['Lab:',handles.ParamSig.Label{iic(jj(1))}]
                ['IC',num2str(iic(jj(1)))],
                ['RT=',num2str(rti(ii(1),jj(1)),3),'ms']
                };
            set(handles.tag_info_point_txt,'string',output_txt)
        end
    end
    
    if get(handles.tag_show_DT,'value')
        dti = handles.MarkersC.dt(:,iic) - repmat(handles.spikes(:),[1 length(iic)]);
        [ii,jj] = find(handles.MarkersC.dt(:,iic)==pos(1)&dti == pos(2));
        
        if ~isempty(ii)
            output_txt= { ['Lab:',handles.ParamSig.Label{iic(jj(1))}]
                ['IC',num2str(iic(jj(1)))],
                ['AT=',num2str(dti(ii(1),jj(1)),3),'ms']
                };
            set(handles.tag_info_point_txt,'string',output_txt)
        end
    end
else
    return
end





function tag_info_point_txt_Callback(hObject, eventdata, handles)
% hObject    handle to tag_info_point_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_info_point_txt as text
%        str2double(get(hObject,'String')) returns contents of tag_info_point_txt as a double



% --- Executes during object creation, after setting all properties.
function tag_info_point_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_info_point_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_dV.
function tag_dV_Callback(hObject, eventdata, handles)
% hObject    handle to tag_dV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_dV


% --- Executes on button press in tag_plot_bipolar.
function tag_plot_bipolar_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_bipolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_bipolar



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


% --- Executes on button press in tag_show_Rw.
function tag_show_Rw_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_Rw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_Rw


% --- Executes on button press in tag_show_RRI.
function tag_show_RRI_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_RRI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_RRI


% --- Executes on button press in tag_hold.
function tag_hold_Callback(hObject, eventdata, handles)
% hObject    handle to tag_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_hold


% --------------------------------------------------------------------
function tag_check_channels_Callback(hObject, eventdata, handles)
% hObject    handle to tag_check_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data.S = handles.S;
data.SNR = handles.SNR;
data.ParamSig = handles.ParamSig;
if isfield(handles,'sig_corr')
    data.sig_corr = handles.sig_corr;
end
if isfield(handles,'Markers')
    data.Markers = handles.Markers;
end
% overwrite with MarkersC
if isfield(handles,'MarkersC')
    data.Markers = handles.MarkersC;
end
if isfield(handles,'MarkersC')
    data.spikes = handles.spikes;
end

% data.handle_ori = handles.hObject_main;
data.handle_ori =hObject;
% Close open windows
aa = findobj(0,'name','Check_Channels_GUI_ATRT');
if ~isempty(aa);
    close(aa);
end

Check_Channels_GUI_ATRT(data);


% --- Executes on button press in tag_show_filter_RT.
function tag_show_filter_RT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_filter_RT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_filter_RT
if handles.tag_show_filter_RT.Value
    handles.tag_show_filter_AT.Value = 0;
    handles.tag_show_filter_RAW.Value = 0;
end
tag_Plot_Callback(hObject, [], handles)



% --- Executes on button press in tag_show_filter_AT.
function tag_show_filter_AT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_filter_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_filter_AT
if handles.tag_show_filter_AT.Value
    handles.tag_show_filter_RT.Value = 0;
    handles.tag_show_filter_RAW.Value = 0;
end
tag_Plot_Callback(hObject, [], handles)

% --- Executes on button press in tag_show_filter_RAW.
function tag_show_filter_RAW_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_filter_RAW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_filter_RAW
if handles.tag_show_filter_RAW.Value
    handles.tag_show_filter_RT.Value = 0;
    handles.tag_show_filter_AT.Value = 0;
end
tag_Plot_Callback(hObject, [], handles)



function tag_HFc_Callback(hObject, eventdata, handles)
% hObject    handle to tag_HFc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_HFc as text
%        str2double(get(hObject,'String')) returns contents of tag_HFc as a double


% --- Executes during object creation, after setting all properties.
function tag_HFc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_HFc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_HFc_AT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_HFc_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_HFc_AT as text
%        str2double(get(hObject,'String')) returns contents of tag_HFc_AT as a double


% --- Executes during object creation, after setting all properties.
function tag_HFc_AT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_HFc_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
