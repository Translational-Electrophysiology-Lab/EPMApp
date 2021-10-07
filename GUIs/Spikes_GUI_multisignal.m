function varargout = Spikes_GUI_multisignal(varargin)
% SPIKES_GUI_MULTISIGNAL MATLAB code for Spikes_GUI_multisignal.fig
%      SPIKES_GUI_MULTISIGNAL, by itself, creates a new SPIKES_GUI_MULTISIGNAL or raises the existing
%      singleton*.
%
%      H = SPIKES_GUI_MULTISIGNAL returns the handle to a new SPIKES_GUI_MULTISIGNAL or the handle to
%      the existing singleton*.
%
%      SPIKES_GUI_MULTISIGNAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPIKES_GUI_MULTISIGNAL.M with the given input arguments.
%
%      SPIKES_GUI_MULTISIGNAL('Property','Value',...) creates a new SPIKES_GUI_MULTISIGNAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Spikes_GUI_multisignal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Spikes_GUI_multisignal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Spikes_GUI_multisignal

% Last Modified by GUIDE v2.5 25-Aug-2021 17:32:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Spikes_GUI_multisignal_OpeningFcn, ...
    'gui_OutputFcn',  @Spikes_GUI_multisignal_OutputFcn, ...
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

% --- Executes just before Spikes_GUI_multisignal is made visible.
function Spikes_GUI_multisignal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Spikes_GUI_multisignal (see VARARGIN)

% Choose default command line output for Spikes_GUI_multisignal
handles.output = hObject;
%% ADD DATA
% data = guidata(varargin{1});
% if isfield(data,'signals')
%     handles.sig = data.signals;
% else
%     handles.sig = data.signals_raw;
% end
% if isfield(data,'spikes')
%     handles.spikes = data.spikes(:).';
% end
% handles.ParamSig = data.ParamSig;
% handles.hObject_main = varargin{1};
% handles.LOG=data.LOG;
% handles.tag_Log=data.tag_Log;
% % handles.tms = [0:size(handles.sig,1)-1]/handles.ParamSig.frequency*1000;
%
% % quality index
% [a,b] = butter(3,[2 30]/(handles.ParamSig.frequency/2));
% xflb = filtfilt(a,b,handles.sig);
% % SNR (estimate)
% qual_ind = 10*log10(sum(xflb.^2)./sum((handles.sig-xflb).^2));
% qual_ind(isnan(qual_ind)) = -inf;
% [~,is] = sort(qual_ind,'descend');
%
%
% List = handles.ParamSig.Label;
% for k = 1:length(List)
%     %     List{k} = ['[',num2str(k),']',List{k}];
%
%     List{k} = ['[',num2str(k),']',List{k},' (',num2str(find(is==k)),')'];
% end
% set(handles.tag_List,'string',List,'value',[1:length(handles.ParamSig.Label)])
% set(handles.tag_samp_ini,'string','1')
% set(handles.tag_samp_fin,'string',length(handles.sig))


% Update handles structure
guidata(hObject, handles);
% UIWAIT makes Spikes_GUI_multisignal wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Spikes_GUI_multisignal_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on selection change in tag_List.
function tag_List_Callback(hObject, eventdata, handles)
% hObject    handle to tag_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_List
cla(handles.tag_ax_sig)
cla(handles.tag_ax_spikes)

tag_Plot_Callback(hObject,[],handles);


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
function tag_Refresh_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hold(handles.tag_ax_spikes,'off')
% hold(handles.tag_ax_sig,'off')
cla(handles.tag_ax_spikes)
cla(handles.tag_ax_sig)
set([handles.tag_ax_sig,handles.tag_ax_spikes],'xlim',[0 1])
% --- Executes on button press in tag_Eliminate.
function tag_Eliminate_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Eliminate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode(handles.figure1);
info_struct = getCursorInfo(dcm_obj);
ic = handles.tag_List.Value;
SpikesUpdated = handles.spikes_all{ic};

ii_out = nan(1,length(info_struct));
spikes_out = nan(1,length(info_struct));
for i = 1:length(info_struct)
    spikes_out(i) = info_struct(i).Position(1);
    switch handles.tag_do_xtick.String{handles.tag_do_xtick.Value}
        case 'Seconds'
            sp = SpikesUpdated/1000;
            spikes_out(i) = spikes_out(i)*1000;
        case 'Minutes'
            sp = SpikesUpdated/1000/60;
            spikes_out(i) = spikes_out(i)*1000*60;
        case 'Hours'
            sp = SpikesUpdated/1000/60/60;
            spikes_out(i) = spikes_out(i)*1000*60*60;
    end
    ii = find(abs(SpikesUpdated-spikes_out(i))<1/handles.ParamSig.frequency);
    if ~isempty(ii)
        ii_out(i)= min(ii); % min only because in theory ii can be a vector
        set(info_struct(i).Target,'visible','off')
    end
end
ii_out(isnan(ii_out))=[];
plot(handles.tag_ax_sig,sp(ii_out),handles.yym,'xr','linewidth',2,'markersize',10)
SpikesUpdated(ii_out)=[];
sp(ii_out)=[];
hold(handles.tag_ax_spikes,'on')
if handles.do_plot_bpm.Value==0
    plot(handles.tag_ax_spikes,sp(2:end),diff(SpikesUpdated),'o-r')
else
    plot(handles.tag_ax_spikes,sp(2:end),60000./diff(SpikesUpdated),'o-r')
end

handles.spikes_all{ic}=SpikesUpdated;
%
set(handles.tag_sdsd,'string',num2str(nanstd(diff(diff(SpikesUpdated))),3));
set(handles.tag_SDNN,'string',num2str(nanstd((diff(SpikesUpdated))),3));


% remove all datatips
delete(findall(gca,'Type','hggroup','HandleVisibility','off'));
% update
guidata(hObject, handles);


% --- Executes on button press in tag_Add.
function tag_Add_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode(handles.figure1);
info_struct = getCursorInfo(dcm_obj);

ic = handles.tag_List.Value;
SpikesUpdated = handles.spikes_all{ic};
spikes_in = nan(1,length(info_struct));
for i = 1:length(info_struct)
    spikes_in(i) = info_struct(i).Position(1);
end

switch handles.tag_do_xtick.String{handles.tag_do_xtick.Value}
    case 'Seconds'
        SpikesUpdated = sort([SpikesUpdated,spikes_in*1000],'ascend');
        sp_plot = SpikesUpdated/1000;
    case 'Minutes'
        SpikesUpdated = sort([SpikesUpdated,spikes_in*1000*60],'ascend');
        sp_plot = SpikesUpdated/1000/60;
    case 'Hours'
        SpikesUpdated = sort([SpikesUpdated,spikes_in*1000*60*60],'ascend');
        sp_plot = SpikesUpdated/1000/60/60;
end

plot(handles.tag_ax_sig,spikes_in,handles.yym,'vg','linewidth',2,'markersize',6,'markerfacecolor','g')
hold(handles.tag_ax_spikes,'on')
if handles.do_plot_bpm.Value==0
    plot(handles.tag_ax_spikes,sp_plot(2:end),diff(SpikesUpdated),'o-r')
else
    plot(handles.tag_ax_spikes,sp_plot(2:end),60000./diff(SpikesUpdated),'o-r')
end
% set(handles.tag_ax_sig,'xtick',sp_plot(2:end),'xticklabel',[],'xgrid','on')
handles.spikes_all{ic}=SpikesUpdated;

set(handles.tag_sdsd,'string',num2str(nanstd(diff(diff(SpikesUpdated))),3));
set(handles.tag_SDNN,'string',num2str(nanstd((diff(SpikesUpdated))),3));


% remove all datatips
delete(findall(gca,'Type','hggroup','HandleVisibility','off'));
% update
guidata(hObject, handles);


% --- Executes on button press in tag_Plot.
function tag_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ic = handles.tag_List.Value;
handles.spikes = handles.spikes_all{ic};
handles.SpikesUpdated = handles.spikes;

xx = get(handles.tag_ax_sig,'xlim');
te = length(handles.sig)/(handles.ParamSig.frequency/1000);
if xx(2)>te
    xx(2)=te;
end



switch handles.tag_do_xtick.String{handles.tag_do_xtick.Value}
    case 'Seconds'
        tms = [0:size(handles.sig,1)-1]/handles.ParamSig.frequency;
        sp_plot = handles.spikes/1000;
        sp_plot_upd = handles.SpikesUpdated/1000;
    case 'Minutes'
        tms = [0:size(handles.sig,1)-1]/handles.ParamSig.frequency/60;
        sp_plot = handles.spikes/1000/60;
        sp_plot_upd = handles.SpikesUpdated/1000/60;
    case 'Hours'
        tms = [0:size(handles.sig,1)-1]/handles.ParamSig.frequency/60/60;
        sp_plot = handles.spikes/1000/60/60;
        sp_plot_upd = handles.SpikesUpdated/1000/60/60;
end
ap2=plot(handles.tag_ax_sig,[1 1]'*sp_plot,[min(handles.sig(:,ic));max(handles.sig(:,ic))]*ones(1,length(handles.spikes)),'--','color',[1 1 1]*.7);
hold(handles.tag_ax_sig,'on')
aa = plot(handles.tag_ax_sig,tms,handles.sig(:,ic),'k');

set([handles.tag_ax_sig,handles.tag_ax_spikes],'xlim',[tms(1) tms(end)])
zoom(handles.tag_ax_sig,'reset')
zoom(handles.tag_ax_spikes,'reset')
zoom(handles.tag_ax_spikes,'out')
zoom(handles.tag_ax_sig,'out')

% vstr = get(handles.tag_List,'string');
% vstr=vstr(ic);
% set(aa(setdiff([1:length(aa)],iiv)),'visible','off')
% ll = legend(aa,vstr);
yy =get(handles.tag_ax_sig,'ylim');
yym = max(.8*max(yy));

if ~isfield(handles,'spikes')
    return
end
plot(handles.tag_ax_sig,sp_plot,yym,'ok','markerfacecolor','r');
clear yy


if handles.do_plot_bpm.Value==0
    y = diff(handles.spikes);
else
    y = 60000./diff(handles.spikes);
end
ap1=plot(handles.tag_ax_spikes,[1 1]'*sp_plot,[min(y)-range(y)*.1;max(y)+range(y)*.1]*ones(1,length(handles.spikes)),'--','color',[1 1 1]*.7);
hold(handles.tag_ax_spikes,'on')
aori = plot(handles.tag_ax_spikes,sp_plot(2:end),y,'o-k','markerfacecolor','r');
set(handles.tag_ax_spikes,'ylim',[min(y)-range(y)*.1 max(y)+range(y)*.1])
if ~isequal(xx,[0 1])
    xlim(handles.tag_ax_sig,xx);
    xlim(handles.tag_ax_spikes,xx);
end
handles.yym = yym;


guidata(hObject, handles);
linkaxes([handles.tag_ax_spikes,handles.tag_ax_sig],'x')

sdsd = nanstd(diff(diff(handles.spikes)));
sdnn = nanstd(diff(handles.spikes));
handles.tag_sdsd.String = num2str(sdsd,3);
handles.tag_SDNN.String = num2str(sdnn,3);


% --- Executes on button press in tag_done.
function tag_done_Callback(hObject, eventdata, handles)
% hObject    handle to tag_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[n,p] = uiputfile(['RestRRV_',handles.tag_List.String{handles.tag_List.Value}]);

signals = handles.sig;
sp_all = handles.spikes_all;
ParamSig = handles.ParamSig;
ic_final = handles.tag_List.Value;
save([p,n],'signals','sp_all','ParamSig','ic_final');



% --- Executes on button press in tag_localize.
function tag_localize_Callback(hObject, eventdata, handles)
% hObject    handle to tag_localize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hold(handles.tag_ax_sig,'off')
plot(handles.tag_ax_sig,nan,nan)
text(.5,.5,' ... Localizing Spikes ....')

hold(handles.tag_ax_spikes,'off')
plot(handles.tag_ax_spikes,nan,nan)
text(.5,.5,' ... Localizing Spikes ....')

set([handles.tag_ax_spikes,handles.tag_ax_sig],'xtick',[],'ytick',[])

ichan = get(handles.tag_List,'value');
%
Hwb = waitbar(0,'Estimating Spikes');
signals = handles.sig(:,ichan);


%% Sinus Rhythm
[spikes] = spikes_detection_SR_CinC(signals,handles.ParamSig.frequency);
% update
handles.spikes_all{ichan}=spikes;
guidata(hObject, handles);
close(Hwb)

tag_Plot_Callback(hObject,[], handles)




function tag_samp_ini_Callback(hObject, eventdata, handles)
% hObject    handle to tag_samp_ini (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_samp_ini as text
%        str2double(get(hObject,'String')) returns contents of tag_samp_ini as a double


% --- Executes during object creation, after setting all properties.
function tag_samp_ini_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_samp_ini (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_samp_fin_Callback(hObject, eventdata, handles)
% hObject    handle to tag_samp_fin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_samp_fin as text
%        str2double(get(hObject,'String')) returns contents of tag_samp_fin as a double


% --- Executes during object creation, after setting all properties.
function tag_samp_fin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_samp_fin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_cut.
function tag_cut_Callback(hObject, eventdata, handles)
% hObject    handle to tag_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ini = str2double(get(handles.tag_samp_ini,'string'));
fin = str2double(get(handles.tag_samp_fin,'string'));
inisamp = round(ini*handles.ParamSig.frequency/1000);
finsamp = round(fin*handles.ParamSig.frequency/1000);

handles.sig = handles.sig(inisamp : finsamp,:);

%
if isfield(handles,'SpikesUpdated')
    handles = rmfield(handles,'SpikesUpdated');
end

if isfield(handles,'spikes')
    handles.spikes(handles.spikes <ini|handles.spikes >fin) = [];
    handles.spikes = handles.spikes-ini+1;
    
    set(handles.tag_sdsd,'string',num2str(nanstd(diff(diff(handles.spikes))),3));
    set(handles.tag_SDNN,'string',num2str(nanstd(diff(handles.spikes)),3));
end
set(handles.tag_samp_ini,'string',num2str(1));
set(handles.tag_samp_fin,'string',num2str(round(size(handles.sig,1)/handles.ParamSig.frequency*1000)));

tms = [0:size(handles.sig,1)-1]/handles.ParamSig.frequency*1000;
hold(handles.tag_ax_sig,'off')
hold(handles.tag_ax_spikes,'off')
plot(handles.tag_ax_sig,tms,handles.sig(:,get(handles.tag_List,'value')));
plot(handles.tag_ax_spikes,nan,nan)
set([handles.tag_ax_spikes,handles.tag_ax_sig],'xlim',[1 size(handles.sig,1)])

set(handles.tag_List,'value',[1:size(handles.sig,2)])
%
set(handles.tag_cut,'backgroundcolor',[.6 .6 .6],'foregroundcolor',[1 0 0])
pause(.5)
set(handles.tag_cut,'backgroundcolor',[1 1 1]*0.941,'foregroundcolor','k')
% update
guidata(hObject,handles)


% --- Executes on button press in tag_do_QRS.
function tag_do_QRS_Callback(hObject, eventdata, handles)
% hObject    handle to tag_do_QRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_do_QRS
if handles.tag_do_QRS.Value == 1
    handles.tag_sensitivity_thr.String = '0.1';
    handles.tag_control.Value = 0;
    handles.tag_VT_check.Value = 0;
    handles.tag_do_QRSAuto.Value = 0;
end

% --- Executes on button press in tag_Eliminate_all.
function tag_Eliminate_all_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Eliminate_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

answ = questdlg('Eliminate all markers currently visualized?');
if isequal(answ,'Yes')
    ic = handles.tag_List.Value;
    
    SpikesUpdated = handles.spikes_all{ic};
    
    xx = get(handles.tag_ax_sig,'xlim');
    switch handles.tag_do_xtick.Value
        case 1
            ii_out = SpikesUpdated/1000<xx(2) & SpikesUpdated/1000>xx(1);
        case 2
            ii_out = SpikesUpdated/1000/60<xx(2) & SpikesUpdated/1000/60>xx(1);
        case 3
            ii_out = SpikesUpdated/1000/60/60<xx(2) & SpikesUpdated/1000/60/60>xx(1);
    end
    
    plot(handles.tag_ax_sig,SpikesUpdated(ii_out),handles.yym,'xr','linewidth',2,'markersize',10)
    SpikesUpdated(ii_out)=[];
    hold(handles.tag_ax_spikes,'on')
    plot(handles.tag_ax_spikes,SpikesUpdated(2:end),diff(SpikesUpdated),'o-r')
    set(handles.tag_ax_sig,'xtick',SpikesUpdated(2:end),'xticklabel',[],'xgrid','on')
    
    handles.spikes_all{ic}=SpikesUpdated;
    
    set(handles.tag_sdsd,'string',num2str(nanstd(diff(diff(handles.spikes))),3));
    set(handles.tag_SDNN,'string',num2str(nanstd((diff(handles.spikes))),3));
    
    % remove all datatips
    delete(findall(gca,'Type','hggroup','HandleVisibility','off'));
    % update
    guidata(hObject, handles);
    
end

function tag_sensitivity_thr_Callback(hObject, eventdata, handles)
% hObject    handle to tag_sensitivity_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_sensitivity_thr as text
%        str2double(get(hObject,'String')) returns contents of tag_sensitivity_thr as a double


% --- Executes during object creation, after setting all properties.
function tag_sensitivity_thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_sensitivity_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_sdsd_Callback(hObject, eventdata, handles)
% hObject    handle to tag_sdsd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_sdsd as text
%        str2double(get(hObject,'String')) returns contents of tag_sdsd as a double


% --- Executes during object creation, after setting all properties.
function tag_sdsd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_sdsd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in tag_cut_view.
function tag_cut_view_Callback(hObject, eventdata, handles)
% hObject    handle to tag_cut_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

xx = get(handles.tag_ax_sig,'xlim');
ini = xx(1);
fin = xx(2);
inisamp = round(ini*handles.ParamSig.frequency/1000);
finsamp = round(fin*handles.ParamSig.frequency/1000);

if inisamp<1
    inisamp = 1;
end
if finsamp>size(handles.sig,1)
    finsamp = size(handles.sig,1);
end


handles.sig = handles.sig(inisamp : finsamp,:);
handles.spikes(handles.spikes <ini|handles.spikes >fin) = [];
handles.spikes = handles.spikes-ini+1;
%
if isfield(handles,'SpikesUpdated')
    handles = rmfield(handles,'SpikesUpdated');
end


sd = diff(handles.spikes);
SD = nanstd(sd)/handles.ParamSig.frequency*1000;
sd(sd>nanmedian(sd)+30 | sd<nanmedian(sd)-30)=[];
SD2 = nanstd(sd)/handles.ParamSig.frequency*1000;
set(handles.tag_sdsd,'string',num2str(SD,3));
if SD2<=1
    set(handles.tag_sdsd,'backgroundcolor','g');
elseif SD2>1&SD2<3
    set(handles.tag_sdsd,'backgroundcolor','y');
else
    set(handles.tag_sdsd,'backgroundcolor','r');
end

set(handles.tag_samp_ini,'string',num2str(1));
set(handles.tag_samp_fin,'string',num2str(round(size(handles.sig,1)/handles.ParamSig.frequency*1000)));

tms = [0:size(handles.sig,1)-1]/handles.ParamSig.frequency*1000;
hold(handles.tag_ax_sig,'off')
hold(handles.tag_ax_spikes,'off')
plot(handles.tag_ax_sig,tms,handles.sig(:,get(handles.tag_List,'value')));
plot(handles.tag_ax_spikes,nan,nan)
set([handles.tag_ax_spikes,handles.tag_ax_sig],'xlim',[1 size(handles.sig,1)])

set(handles.tag_List,'value',[1:size(handles.sig,2)])
%
set(handles.tag_cut,'backgroundcolor',[.6 .6 .6],'foregroundcolor',[1 0 0])
pause(.5)
set(handles.tag_cut,'backgroundcolor',[1 1 1]*0.941,'foregroundcolor','k')
% update
guidata(hObject,handles)


% --- Executes on button press in tag_localize_in_interval.
function tag_localize_in_interval_Callback(hObject, eventdata, handles)
% hObject    handle to tag_localize_in_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Same as Localize, but within an interval, keeping the spikes previously estimated

tt = get(handles.tag_ax_sig,'xlim');
if tt(1)==0
    helpdlg('Please select a specific interval')
    return
end

if tt(1)<1, tt(1)=1; end
if tt(end)>size(handles.sig,1)/handles.ParamSig.frequency*1000, tt(end)=size(handles.sig,1)/handles.ParamSig.frequency*1000; end
if isfield(handles,'spikes');
    ii = handles.spikes<=tt(2)&handles.spikes>=tt(1);
    handles.spikes(ii)=[];
end
clear ii
ttsamp = round(tt/1000*handles.ParamSig.frequency);
% ttsamp(ttsamp<1 | ttsamp>size(handles.sig,1))=[];
hold(handles.tag_ax_sig,'off')
plot(handles.tag_ax_sig,nan,nan)
text(.5,.5,' ... Localizing Spikes ....')

hold(handles.tag_ax_spikes,'off')
plot(handles.tag_ax_spikes,nan,nan)
text(.5,.5,' ... Localizing Spikes ....')

set([handles.tag_ax_spikes,handles.tag_ax_sig],'xtick',[],'ytick',[])

ichan = get(handles.tag_List,'value');
handles.LOG = [handles.LOG,['- Finding spikes using ',num2str(length(ichan)),' channels']];
set(handles.tag_Log,'string',handles.LOG)
%
Hwb = waitbar(0,'Estimating Spikes');

%% only interval visualized
signals = handles.sig(ttsamp(1):ttsamp(2),ichan);
if get(handles.tag_control,'value')
    do_control = get(handles.tag_control,'value');
    sp_all = cell(1,size(signals,2));
    sp_all_D = cell(1,size(signals,2));
    sensthr = str2num(get(handles.tag_sensitivity_thr,'string'));
    for ind = 1:size(signals,2)
        waitbar(ind/size(signals,2),Hwb);
        if isequal(handles.ParamSig.geoname,'MEA')
            minsp = 20;
        else
            minsp = 150;
        end
        [sp] = find_pacing_spikes_mo(signals(:,ind),handles.ParamSig.frequency,sensthr,do_control,minsp);
        sp_all{ind} = sp;
        sp_all_D{ind} = diff(sp);
    end
    close(Hwb)
    
    
    [SD,iiSD] = min(cellfun(@nanstd,sp_all_D));
    spikes = sp_all{iiSD};
    % Update spikes
    spikes = (spikes)/handles.ParamSig.frequency*1000 + tt(1)-1;
    spikes = unique([handles.spikes spikes]);
    SD =  SD/handles.ParamSig.frequency*1000;
    set(handles.tag_List,'value',ichan(iiSD));
    %% plot
    hold(handles.tag_ax_spikes,'off')
    hold(handles.tag_ax_sig,'off')
    tms = [0:size(handles.sig,1)-1]/handles.ParamSig.frequency*1000;
    %     plot(handles.tag_ax_sig,tms,signals(:,iiSD)),
    plot(handles.tag_ax_sig,tms,handles.sig(:,ichan(iiSD))),
    hold(handles.tag_ax_sig,'on')
    
    hh = get(handles.tag_ax_sig,'ylim');
    for i = 1:length(spikes)
        plot(handles.tag_ax_sig,[1,1]*spikes(i),hh,'--k')
    end
    
    
    plot(handles.tag_ax_spikes,spikes(2:end),diff(spikes),'o-')
    %     hold(handles.tag_ax_sig,'on')
    %     yym = max(max(signals(round(spikes*handles.ParamSig.frequency/1000),iiSD)));
    yym = max(max(signals(:,iiSD)));
    
    
    plot(handles.tag_ax_sig,spikes,yym,'ob');
    
    %     set([handles.tag_ax_sig],'xlim',[tms(1) tms(end)]);
    linkaxes([handles.tag_ax_sig,handles.tag_ax_spikes],'x')
    set([handles.tag_ax_sig,handles.tag_ax_spikes],'xtick',spikes)
    set(handles.tag_ax_spikes,'xticklabel',[1:length(spikes)],'ygrid','on')
    set(handles.tag_ax_sig,'xticklabel','')
    title(handles.tag_ax_sig,['Median CL = ',num2str(nanmedian(diff(spikes)),4),'ms - HeartBeats #',num2str(length(spikes))])
    ylabel(handles.tag_ax_spikes,'[ms]')
    
    yy =get(handles.tag_ax_sig,'ylim');
    handles.yym = yym;
    
    sd = sp_all_D{iiSD};
    sd(sd>nanmedian(sd)+30 | sd<nanmedian(sd)-30)=[];
    SD2 = nanstd(sd)/handles.ParamSig.frequency*1000;
    set(handles.tag_sdsd,'string',num2str(SD,3));
    
    if SD2<=1
        set(handles.tag_sdsd,'backgroundcolor','g');
    elseif SD2>1&SD2<3
        set(handles.tag_sdsd,'backgroundcolor','y');
    else
        set(handles.tag_sdsd,'backgroundcolor','r');
    end
    
    
    % update
    handles.spikes=spikes;
    % end
    handles.SpikesUpdated=spikes;
    guidata(hObject, handles);
    
    %% Sinus Rhythm
elseif get(handles.tag_do_QRS,'value')
    
    sp_all = cell(1,size(signals,2));
    sp_all_D = cell(1,size(signals,2));
    x = signals;
    
    waitbar(1/3,Hwb);
    [a,b] = butter(3,[2 30]/(handles.ParamSig.frequency/2));
    xflb = filtfilt(a,b,x);
    
    Wmed = 25;%50
    xf = medfilt1(diff(xflb).^2,round(Wmed/1000*handles.ParamSig.frequency));
    %     xf = medfilt1(diff(x).^2,50);
    
    sensthr = str2num(get(handles.tag_sensitivity_thr,'string'));
    h = xf>repmat(max(xf)*sensthr,[size(xf,1) 1]);
    h = imclose(h,ones(round(Wmed/1000*handles.ParamSig.frequency),size(h,2)));
    %     L = 10;
    %     xd = diff(filtfilt(hanning(L),sum(hanning(L)),x));
    %     xd = diff(xflb);
    xd = diff(x);
    
    xd2 = diff(xd);
    for k = 1:size(h,2)
        waitbar(1/3 + 2/3*k/size(signals,2),Hwb);
        ih = find(diff(h(:,k))>.99);
        eh = find(diff(h(:,k))<-.99); % modify to find VT/VF 02/02/2015
        while eh(1)<ih(1)
            eh(1)=[];
        end
        Lm = min(length(ih),length(eh));
        ih = ih(1:Lm); eh = eh(1:Lm);
        clear Lm
        % -
        sp = nan(1,length(ih));
        for j = 1:length(ih)
            H = (ih(j)-round(40/1000*handles.ParamSig.frequency)):(eh(j)+round(40/1000*handles.ParamSig.frequency)); %
            %             H = (ih(j)-round(20/1000*handles.ParamSig.frequency)):(eh(j)+round(20/1000*handles.ParamSig.frequency)); % modify to find VT/VF 02/02/2015
            H(H<1 | H>size(xd2,1))=[];
            %             ii = find([(xd(H(1:length(H)-1),k).*xd(H(2:length(H)),k))<=0 & xd2(H(1:end-1),k)<0])+H(1)+1;
            ii = find((xd(H(1:length(H)-1),k).*xd(H(2:length(H)),k))<=0) + 1;
            if length(ii)>1
                %                 [~,id] = min(xd2(ii,k)); % min of second derivative: it
                %                 works for clean and "spiky" QRS
                %                 x2 = detrend(x(H,k),'constant');
                [~,id] = max(abs(x(H(ii)))); %
                ii = ii(id);
            end
            if ~isempty(ii)
                sp_all{k}(j) = H(ii);
            end
            
            %             figure,plot(H,xf(H,k)),hold on,plot(H,x(H,k),'r'),plot(sp_all{k}(j),x(sp_all{k}(j),k),'or')
        end
        sp_all{k}(sp_all{k}==0 | isnan(sp_all{k}) | [false diff(sp_all{k})<1]) = [];
        sp_all_D{k} = diff(sp_all{k});
        
        %     figure,plot(xf(:,k)),hold on,plot(x(:,k),'r'),plot(sp_all{k}(:),x(sp_all{k}(:),k),'or')
    end
    
    close(Hwb)
    
    %     iiok = find(~cellfun(@isempty,sp_all_D)&cellfun(@length,sp_all_D)>5);
    iiok = find(~cellfun(@isempty,sp_all_D));
    
    if ~isempty(iiok)
        [SD,iiSD] = min(cellfun(@nanstd,sp_all_D(iiok)));
        iiSD = iiok(iiSD);
        spikes = sp_all{iiSD};
        %     spikes = spikes/handles.ParamSig.frequency*1000;
        % Update spikes
        spikes = (spikes)/handles.ParamSig.frequency*1000 + tt(1)-1;
        if isfield(handles,'spikes')
            spikes = unique([handles.spikes spikes]);
        end
        SD =  SD/handles.ParamSig.frequency*1000;
        set(handles.tag_List,'value',ichan(iiSD));
        %% plot
        tms = [1:size(handles.sig(:,ichan(iiSD)),1)]/handles.ParamSig.frequency*1000;
        hold(handles.tag_ax_sig,'off')
        %     plot(handles.tag_ax_sig,tms,signals(:,iiSD)),
        plot(handles.tag_ax_sig,tms,handles.sig(:,ichan(iiSD))),
        hold(handles.tag_ax_sig,'on')
        
        hh = get(handles.tag_ax_sig,'ylim');
        for i = 1:length(spikes)
            plot(handles.tag_ax_sig,[1,1]*spikes(i),hh,'--k')
        end
        
        plot(handles.tag_ax_spikes,spikes(2:end),diff(spikes),'o-')
        hold(handles.tag_ax_sig,'on')
        yym = get(handles.tag_ax_sig,'ylim');
        plot(handles.tag_ax_sig,spikes,yym(2),'ob');
        
        set([handles.tag_ax_sig],'xlim',[tms(1) tms(end)]);
        linkaxes([handles.tag_ax_sig,handles.tag_ax_spikes],'x')
        set([handles.tag_ax_sig,handles.tag_ax_spikes],'xtick',spikes)
        set(handles.tag_ax_spikes,'xticklabel',[1:length(spikes)],'ygrid','on')
        set(handles.tag_ax_sig,'xticklabel','')
        title(handles.tag_ax_sig,['Median CL = ',num2str(nanmedian(diff(spikes)),4),'ms - HeartBeats #',num2str(length(spikes))])
        ylabel(handles.tag_ax_spikes,'[ms]')
        
        yy =get(handles.tag_ax_sig,'ylim');
        handles.yym = yym;
        
        sd = sp_all_D{iiSD};
        sd(sd>nanmedian(sd)+30 | sd<nanmedian(sd)-30)=[];
        SD2 = nanstd(sd)/handles.ParamSig.frequency*1000;
        set(handles.tag_sdsd,'string',num2str(SD,3));
        
        if SD2<=1
            set(handles.tag_sdsd,'backgroundcolor','g');
        elseif SD2>1&SD2<3
            set(handles.tag_sdsd,'backgroundcolor','y');
        else
            set(handles.tag_sdsd,'backgroundcolor','r');
        end
        %% Correction
        %     [spikes] = spikes_correction_SR(signals,spikes,handles.ParamSig,1);
        
        % update
        handles.spikes=spikes;
        handles.SpikesUpdated=spikes;
        guidata(hObject, handles);
        
    else
        helpdlg('No spike detected')
        pause(1)
        close
    end
    
    
elseif get(handles.tag_VT_check,'value')
    
    sp_all = cell(1,size(signals,2));
    sp_all_D = cell(1,size(signals,2));
    x = signals;
    
    % waitbar(1/3,Hwb);
    [a,b] = butter(3,[2 30]/(handles.ParamSig.frequency/2));
    xflb = filtfilt(a,b,x);
    
    Wmed = 25;%50
    xf = medfilt1(diff(xflb),round(Wmed/1000*handles.ParamSig.frequency));
    
    sensthr = str2num(get(handles.tag_sensitivity_thr,'string'));
    h = (xf(1:end-1).*xf(2:end))<0 & abs(x(3:end))>sensthr*prctile(abs(x),98);
    
    for k = 1:size(h,2)
        
        ih = find(h);
        if mean(xflb(h)>0)>.5
            ih(xflb(ih)<0) = [];
        else
            ih(xflb(ih)>0) = [];
        end
        
        sp_all{k} = ih(:).';
        sp_all{k}(sp_all{k}==0 | isnan(sp_all{k}) | [false diff(sp_all{k})<1]) = [];
        sp_all_D{k} = diff(sp_all{k});
    end
    
    close(Hwb)
    iiok = find(~cellfun(@isempty,sp_all_D));
    
    if ~isempty(iiok)
        [SD,iiSD] = min(cellfun(@nanstd,sp_all_D(iiok)));
        iiSD = iiok(iiSD);
        spikes = sp_all{iiSD};
        %     spikes = spikes/handles.ParamSig.frequency*1000;
        % Update spikes
        spikes = (spikes)/handles.ParamSig.frequency*1000 + tt(1)-1;
        spikes = unique([handles.spikes spikes]);
        
        SD =  SD/handles.ParamSig.frequency*1000;
        set(handles.tag_List,'value',ichan(iiSD));
        %% plot
        tms = [1:size(handles.sig(:,ichan(iiSD)),1)]/handles.ParamSig.frequency*1000;
        hold(handles.tag_ax_sig,'off')
        %     plot(handles.tag_ax_sig,tms,signals(:,iiSD)),
        plot(handles.tag_ax_sig,tms,handles.sig(:,ichan(iiSD))),
        hold(handles.tag_ax_sig,'on')
        
        hh = get(handles.tag_ax_sig,'ylim');
        for i = 1:length(spikes)
            plot(handles.tag_ax_sig,[1,1]*spikes(i),hh,'--k')
        end
        
        plot(handles.tag_ax_spikes,spikes(2:end),diff(spikes),'o-')
        hold(handles.tag_ax_sig,'on')
        yym = get(handles.tag_ax_sig,'ylim');
        plot(handles.tag_ax_sig,spikes,yym(2),'ob');
        
        set([handles.tag_ax_sig],'xlim',[tms(1) tms(end)]);
        linkaxes([handles.tag_ax_sig,handles.tag_ax_spikes],'x')
        set([handles.tag_ax_sig,handles.tag_ax_spikes],'xtick',spikes)
        set(handles.tag_ax_spikes,'xticklabel',[1:length(spikes)],'ygrid','on')
        set(handles.tag_ax_sig,'xticklabel','')
        title(handles.tag_ax_sig,['Median CL = ',num2str(nanmedian(diff(spikes)),4),'ms - HeartBeats #',num2str(length(spikes))])
        ylabel(handles.tag_ax_spikes,'[ms]')
        
        yy =get(handles.tag_ax_sig,'ylim');
        handles.yym = yym;
        
        sd = sp_all_D{iiSD};
        sd(sd>nanmedian(sd)+30 | sd<nanmedian(sd)-30)=[];
        SD2 = nanstd(sd)/handles.ParamSig.frequency*1000;
        set(handles.tag_sdsd,'string',num2str(SD,3));
        
        if SD2<=1
            set(handles.tag_sdsd,'backgroundcolor','g');
        elseif SD2>1&SD2<3
            set(handles.tag_sdsd,'backgroundcolor','y');
        else
            set(handles.tag_sdsd,'backgroundcolor','r');
        end
        %% Correction
        %     [spikes] = spikes_correction_SR(signals,spikes,handles.ParamSig,1);
        
        % update
        handles.spikes=spikes;
        handles.SpikesUpdated=spikes;
        guidata(hObject, handles);
        
    else
        helpdlg('No spike detected')
        pause(1)
        close
    end
    
    
    
end

set([handles.tag_ax_sig,handles.tag_ax_spikes],'xlim',[0 size(handles.sig,1)/handles.ParamSig.frequency*1000])
zoom(handles.tag_ax_sig,'reset')
zoom(handles.tag_ax_spikes,'reset')
zoom(handles.tag_ax_spikes,'out')
zoom(handles.tag_ax_sig,'out')
set([handles.tag_ax_sig],'xlim',[tt(1) tt(end)]);

guidata(hObject, handles);


% --- Executes on button press in tag_VT_check.
function tag_VT_check_Callback(hObject, eventdata, handles)
% hObject    handle to tag_VT_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_VT_check
if handles.tag_VT_check.Value == 1
    handles.tag_sensitivity_thr.String = '0.5';
    handles.tag_do_QRS.Value = 0;
    handles.tag_do_QRSAuto.Value = 0;
    handles.tag_control.Value = 0;
end


% --- Executes on button press in tag_do_QRSAuto.
function tag_do_QRSAuto_Callback(hObject, eventdata, handles)
% hObject    handle to tag_do_QRSAuto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.tag_do_QRSAuto.Value == 1
    handles.tag_sensitivity_thr.String = 'NA';
    handles.tag_do_QRS.Value = 0;
    handles.tag_control.Value = 0;
    handles.tag_VT_check.Value = 0;
end

% --- Executes on button press in tag_control.
function tag_control_Callback(hObject, eventdata, handles)
% hObject    handle to tag_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_control
if handles.tag_control.Value == 1
    handles.tag_sensitivity_thr.String = '0.7';
    handles.tag_do_QRS.Value = 0;
    handles.tag_VT_check.Value = 0;
    handles.tag_do_QRSAuto.Value = 0;
end


% --- Executes on button press in do_plot_bpm.
function do_plot_bpm_Callback(hObject, eventdata, handles)
% hObject    handle to do_plot_bpm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of do_plot_bpm


% --- Executes on selection change in tag_do_xtick.
function tag_do_xtick_Callback(hObject, eventdata, handles)
% hObject    handle to tag_do_xtick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_do_xtick contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_do_xtick


% --- Executes during object creation, after setting all properties.
function tag_do_xtick_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_do_xtick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_load.
function tag_load_Callback(hObject, eventdata, handles)
% hObject    handle to tag_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f,p] = uigetfile('*.mat');
V = load([p,f]);
handles.sig = V.signals;
handles.spikes_all = V.sp_all;
handles.ParamSig = V.ParamSig;
% handles.ParamSig.Label = V.pt90;
% handles.RRV = V.RestRRV;
set(handles.tag_List,'string',handles.ParamSig.Label,'value',[1:length(handles.ParamSig.Label)]);

guidata(hObject, handles);



function tag_SDNN_Callback(hObject, eventdata, handles)
% hObject    handle to tag_SDNN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_SDNN as text
%        str2double(get(hObject,'String')) returns contents of tag_SDNN as a double


% --- Executes during object creation, after setting all properties.
function tag_SDNN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_SDNN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
