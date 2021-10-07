function varargout = Spikes_GUI_standalone(varargin)

%% Spikes_GUI_standalone(signals(:,iis(1:20)),ParamSig.frequency,spikes,ParamSig.Label(iis(1:20)),do_wait);

% SPIKES_GUI_STANDALONE MATLAB code for Spikes_GUI_standalone.fig
%      SPIKES_GUI_STANDALONE, by itself, creates a new SPIKES_GUI_STANDALONE or raises the existing
%      singleton*.
%
%      H = SPIKES_GUI_STANDALONE returns the handle to a new SPIKES_GUI_STANDALONE or the handle to
%      the existing singleton*.
%
%      SPIKES_GUI_STANDALONE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPIKES_GUI_STANDALONE.M with the given input arguments.
%
%      SPIKES_GUI_STANDALONE('Property','Value',...) creates a new SPIKES_GUI_STANDALONE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Spikes_GUI_standalone_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Spikes_GUI_standalone_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Spikes_GUI_standalone

% Last Modified by GUIDE v2.5 24-Jan-2016 18:46:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Spikes_GUI_standalone_OpeningFcn, ...
    'gui_OutputFcn',  @Spikes_GUI_standalone_OutputFcn, ...
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

% --- Executes just before Spikes_GUI_standalone is made visible.
function Spikes_GUI_standalone_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Spikes_GUI_standalone (see VARARGIN)

% Choose default command line output for Spikes_GUI_standalone
handles.output = hObject;
%% ADD DATA
handles.sig = varargin{1};
if length(varargin)==5
    do_wait = varargin{5};
    handles.ParamSig.Label = varargin{4};
    handles.spikes = varargin{3};
    handles.ParamSig.frequency = varargin{2};
end

if length(varargin)==4
    do_wait = 0;
    handles.ParamSig.Label = varargin{4};
    handles.spikes = varargin{3};
    handles.ParamSig.frequency = varargin{2};
end

if length(varargin)==3
    do_wait = 0;
    handles.ParamSig.Label = cell(size(handles.sig,2),1);
    handles.spikes = varargin{3};
    handles.ParamSig.frequency = varargin{2};
end

if length(varargin)==2
    do_wait = 0;
    handles.ParamSig.frequency = varargin{2};
    handles.ParamSig.Label = cell(size(handles.sig,2),1);
    for i = 1:size(handles.sig,2)
        handles.ParamSig.Label(i) = {['IC-',num2str(i)]};
    end
    
end

if length(varargin)==1
    A = inputdlg('Introduce Sampling Rate in Hz');
    handles.ParamSig.frequency = str2double(A);clear A
    handles.ParamSig.Label = cell(size(handles.sig,2),1);
    for i = 1:size(handles.sig,2)
        handles.ParamSig.Label(i) = {['IC-',num2str(i)]};
    end
    do_wait = 0;
end

if isempty(handles.ParamSig.frequency)
    A = inputdlg('Introduce Sampling Rate in Hz');
    handles.ParamSig.frequency = str2double(A);clear A
end
if isempty(handles.ParamSig.Label)
    handles.ParamSig.Label = cell(size(handles.sig,2),1);
    for i = 1:size(handles.sig,2)
        handles.ParamSig.Label(i) = {['IC-',num2str(i)]};
    end
end
%%
% quality index
[a,b] = butter(3,[2 30]/(handles.ParamSig.frequency/2));
xflb = filtfilt(a,b,handles.sig);
% SNR (estimate)
qual_ind = 10*log10(sum(xflb.^2)./sum((handles.sig-xflb).^2));
qual_ind(isnan(qual_ind)) = -inf;
[~,is] = sort(qual_ind,'descend');


List = handles.ParamSig.Label;
for k = 1:length(List)
    List{k} = ['[',num2str(k),']',List{k},' (',num2str(find(is==k)),')'];
end
set(handles.tag_List,'string',List,'value',[1:length(handles.ParamSig.Label)])
set(handles.tag_samp_ini,'string','1')
set(handles.tag_samp_fin,'string',length(handles.sig))


handles.do_wait = do_wait;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes Spikes_GUI_standalone wait for user response (see UIRESUME)
if do_wait
    uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = Spikes_GUI_standalone_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = [];

% --- Executes on selection change in tag_List.
function tag_List_Callback(hObject, eoneventdata, handles)
% hObject    handle to tag_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_List


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
hold(handles.tag_ax_spikes,'off')
hold(handles.tag_ax_sig,'off')
% tms = [0:size(handles.sig,1)-1]/handles.ParamSig.frequency*1000;
% xlim([tms(1) tms(end)])

% --- Executes on button press in tag_Eliminate.
function tag_Eliminate_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Eliminate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode(handles.figure1);
info_struct = getCursorInfo(dcm_obj);

if isfield(handles,'SpikesUpdated');
    SpikesUpdated = handles.SpikesUpdated;
else
    SpikesUpdated = handles.spikes;
end
ii_out = nan(1,length(info_struct));
for i = 1:length(info_struct)
    spikes_out(i) = info_struct(i).Position(1);
    ii = find(SpikesUpdated==spikes_out(i));
    if ~isempty(ii)
        ii_out(i)= min(ii); % min only because in theory ii can be a vector
        set(info_struct(i).Target,'visible','off')
    end
end
ii_out(isnan(ii_out))=[];
plot(handles.tag_ax_sig,SpikesUpdated(ii_out),handles.yym,'xr','linewidth',2,'markersize',10)
SpikesUpdated(ii_out)=[];
hold(handles.tag_ax_spikes,'on')
plot(handles.tag_ax_spikes,SpikesUpdated(2:end),diff(SpikesUpdated),'o-r')
set(handles.tag_ax_sig,'xtick',SpikesUpdated(2:end),'xticklabel',[],'xgrid','on')

handles.SpikesUpdated=SpikesUpdated;
%
sd = diff(SpikesUpdated);
SD = nanstd(sd)/handles.ParamSig.frequency*1000;
sd(sd>nanmedian(sd)+30 | sd<nanmedian(sd)-30)=[];
SD2 = nanstd(sd)/handles.ParamSig.frequency*1000;
set(handles.tag_std,'string',num2str(SD,3));
if SD2<=1
    set(handles.tag_std,'backgroundcolor','g');
elseif SD2>1&SD2<3
    set(handles.tag_std,'backgroundcolor','y');
else
    set(handles.tag_std,'backgroundcolor','r');
end

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

if isfield(handles,'SpikesUpdated');
    SpikesUpdated = handles.SpikesUpdated;
else
    SpikesUpdated = handles.spikes;
end

for i = 1:length(info_struct)
    spikes_in(i) = info_struct(i).Position(1);
end
SpikesUpdated = sort([SpikesUpdated,spikes_in],'ascend');
SpikesUpdated(isnan(SpikesUpdated)) = [];
plot(handles.tag_ax_sig,spikes_in,handles.yym,'vg','linewidth',2,'markersize',6,'markerfacecolor','g')
hold(handles.tag_ax_spikes,'on')
plot(handles.tag_ax_spikes,SpikesUpdated(2:end),diff(SpikesUpdated),'o-r')
set(handles.tag_ax_sig,'xtick',SpikesUpdated(2:end),'xticklabel',[],'xgrid','on')
handles.SpikesUpdated=SpikesUpdated;

sd = diff(SpikesUpdated);
SD = nanstd(sd)/handles.ParamSig.frequency*1000;
sd(sd>nanmedian(sd)+30 | sd<nanmedian(sd)-30)=[];
SD2 = nanstd(sd)/handles.ParamSig.frequency*1000;
set(handles.tag_std,'string',num2str(SD,3));
if SD2<=1
    set(handles.tag_std,'backgroundcolor','g');
elseif SD2>1&SD2<3
    set(handles.tag_std,'backgroundcolor','y');
else
    set(handles.tag_std,'backgroundcolor','r');
end
% remove all datatips
delete(findall(gca,'Type','hggroup','HandleVisibility','off'));
% update
guidata(hObject, handles);


% --- Executes on button press in tag_Plot.
function tag_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

xx = get(handles.tag_ax_sig,'xlim');
te = length(handles.sig)/(handles.ParamSig.frequency/1000);
if xx(2)>te
    xx(2)=te;
end

tms = [0:size(handles.sig,1)-1]/handles.ParamSig.frequency*1000;
aa = plot(handles.tag_ax_sig,tms,handles.sig);
set([handles.tag_ax_sig,handles.tag_ax_spikes],'xlim',[tms(1) tms(end)])
zoom(handles.tag_ax_sig,'reset')
zoom(handles.tag_ax_spikes,'reset')
zoom(handles.tag_ax_spikes,'out')
zoom(handles.tag_ax_sig,'out')

iiv = get(handles.tag_List,'value');
vstr = get(handles.tag_List,'string');
vstr=vstr(iiv);
set(aa(setdiff([1:length(aa)],iiv)),'visible','off')
ll = legend(aa(iiv),vstr);
hold(handles.tag_ax_sig,'on')
yy =get(handles.tag_ax_sig,'ylim');
yym = max(.8*max(yy));

if ~isfield(handles,'spikes')
    return
end
plot(handles.tag_ax_sig,handles.spikes,yym,'ob');
clear yy

aori = plot(handles.tag_ax_spikes,handles.spikes(2:end),diff(handles.spikes),'.-');
% linkaxes([handles.tag_ax_spikes,handles.tag_ax_sig],'x')
if ~isequal(xx,[0 1])
    %     set(handles.tag_ax_sig,'xlim',xx);
    xlim(handles.tag_ax_sig,xx);
    xlim(handles.tag_ax_spikes,xx);
end
if isfield(handles,'SpikesUpdated')
    hold(handles.tag_ax_spikes,'on')
    if ~isempty(setdiff(handles.spikes,handles.SpikesUpdated))
        plot(handles.tag_ax_sig,setdiff(handles.spikes,handles.SpikesUpdated),yym,'xr','linewidth',2,'markersize',10)
    end
    if ~isempty(setdiff(handles.SpikesUpdated,handles.spikes))
        plot(handles.tag_ax_sig,setdiff(handles.SpikesUpdated,handles.spikes),yym,'vg','linewidth',2,'markersize',6,'markerfacecolor','g')
    end
    plot(handles.tag_ax_spikes,handles.SpikesUpdated(2:end),diff(handles.SpikesUpdated),'o-r')
    set(handles.tag_ax_sig,'xtick',handles.SpikesUpdated(2:end),'xticklabel',[],'xgrid','on')
    set(handles.tag_ax_spikes,'xtick',handles.SpikesUpdated(2:end),'xticklabel',[1:length(handles.SpikesUpdated(2:end))],'ygrid','on')
    
    set(aori,'color',[.8 .8 .8])
else
    set(handles.tag_ax_spikes,'xtick',handles.spikes(2:end),'xticklabel',[1:length(handles.spikes(2:end))],'ygrid','on')
    set(handles.tag_ax_sig,'xtick',handles.spikes(2:end),'xticklabel',[],'xgrid','on')
end

handles.yym = yym;
% update
guidata(hObject, handles);
linkaxes([handles.tag_ax_spikes,handles.tag_ax_sig],'x')


% --- Executes on button press in tag_done.
function tag_done_Callback(hObject, eventdata, handles)
% hObject    handle to tag_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'SpikesUpdated')
    data.spikes = handles.SpikesUpdated;
else
    data.spikes = handles.spikes;
end
assignin('base','DataFromGui',data);
% % update
guidata(hObject,handles);
close(handles.figure1)

% if ~handles.do_wait
%     close(handles.figure1)
% end

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

if get(handles.tag_do_QRS,'value')==0
    
    do_control = get(handles.tag_control,'value');
    sp_all = cell(1,size(signals,2));
    sp_all_D = cell(1,size(signals,2));
    sensthr = str2num(get(handles.tag_sensitivity_thr,'string'));
    for ind = 1:size(signals,2)
        waitbar(ind/size(signals,2),Hwb);
        
        minsp = 150;
        
        [sp] = find_pacing_spikes_mo(signals(:,ind),handles.ParamSig.frequency,sensthr,do_control,minsp);
        sp_all{ind} = sp;
        sp_all_D{ind} = diff(sp);
        %         sig3_all(:,ind) = butterworthfilter(sig2,handles.ParamSig.frequency,[0.5 300]);
    end
    close(Hwb)
    
    
    [SD,iiSD] = min(cellfun(@nanstd,sp_all_D));
    spikes = sp_all{iiSD};
    spikes = spikes/handles.ParamSig.frequency*1000;
    SD =  SD/handles.ParamSig.frequency*1000;
    set(handles.tag_List,'value',ichan(iiSD));
    %% plot
    hold(handles.tag_ax_spikes,'off')
    hold(handles.tag_ax_sig,'off')
    tms = [0:size(signals,1)-1]/handles.ParamSig.frequency*1000;
    plot(handles.tag_ax_sig,tms,signals(:,iiSD)),
    hold(handles.tag_ax_sig,'on')
    
    hh = get(handles.tag_ax_sig,'ylim');
    for i = 1:length(spikes)
        plot(handles.tag_ax_sig,[1,1]*spikes(i),hh,'--k')
    end
    
    
    plot(handles.tag_ax_spikes,spikes(2:end),diff(spikes),'o-')
    %     hold(handles.tag_ax_sig,'on')
    %     yym = max(max(signals(round(spikes*handles.ParamSig.frequency/1000),iiSD)));
    yym = max(max(signals(:,iiSD)));
    
    if isempty(spikes)
        helpdlg('No spike detected')
        pause(1)
        close
    else
        plot(handles.tag_ax_sig,spikes,yym,'ob');
    end
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
    set(handles.tag_std,'string',num2str(SD,3));
    
    if SD2<=1
        set(handles.tag_std,'backgroundcolor','g');
    elseif SD2>1&SD2<3
        set(handles.tag_std,'backgroundcolor','y');
    else
        set(handles.tag_std,'backgroundcolor','r');
    end
    
    
    % update
    % if ~isfield(handles,'spikes')
    handles.spikes=spikes;
    % end
    handles.SpikesUpdated=spikes;
    guidata(hObject, handles);
    
    %% Sinus Rhythm
else
    
    sp_all = cell(1,size(signals,2));
    sp_all_D = cell(1,size(signals,2));
    x = signals;
    
    waitbar(1/3,Hwb);
    % mo 05/2015
    %     [a,b] = cheby1(7,0.5,35/(handles.ParamSig.frequency/2));
    [a,b] = butter(3,[2 30]/(handles.ParamSig.frequency/2));
    xflb = filtfilt(a,b,x);
    
    
    xf = medfilt1(diff(xflb).^2,round(50/1000*handles.ParamSig.frequency));
    %     xf = medfilt1(diff(x).^2,50);
    
    sensthr = str2num(get(handles.tag_sensitivity_thr,'string'));
    h = xf>repmat(max(xf)*sensthr,[size(xf,1) 1]);
    h = imclose(h,ones(round(50/1000*handles.ParamSig.frequency),size(h,2)));
    L = 10;
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
            %             H = ih(j) + [-20:40];
            H = (ih(j)-round(20/1000*handles.ParamSig.frequency)):(eh(j)+round(20/1000*handles.ParamSig.frequency)); % modify to find VT/VF 02/02/2015
            H(H<1 | H>size(xd2,1))=[];
            ii= find([(xd(H(1:length(H)-1),k).*xd(H(2:length(H)),k))<=0 & xd2(H(1:end-1),k)<0])+H(1)+1;
            if length(ii)>1
                %                 [~,id] = min(xd2(ii,k)); % min of second derivative: it
                %                 works for clean and "spiky" QRS
                x2 = detrend(x(H,k));
                [~,id] = max(abs(x2(ii-H(1)))); %
                ii = ii(id);
            end
            if ~isempty(ii)
                sp_all{k}(j) = ii-1;
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
        spikes = spikes/handles.ParamSig.frequency*1000;
        SD =  SD/handles.ParamSig.frequency*1000;
        set(handles.tag_List,'value',ichan(iiSD));
        %% plot
        tms = [1:size(signals,1)]/handles.ParamSig.frequency*1000;
        hold(handles.tag_ax_sig,'off')
        plot(handles.tag_ax_sig,tms,signals(:,iiSD)),
        hold(handles.tag_ax_sig,'on')
        
        hh = get(handles.tag_ax_sig,'ylim');
        for i = 1:length(spikes)
            plot(handles.tag_ax_sig,[1,1]*spikes(i),hh,'--k')
        end
        
        plot(handles.tag_ax_spikes,spikes(2:end),diff(spikes),'o-')
        hold(handles.tag_ax_sig,'on')
        yym = max(max(signals(round(spikes*handles.ParamSig.frequency/1000),:)));
        plot(handles.tag_ax_sig,spikes,yym,'ob');
        
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
        set(handles.tag_std,'string',num2str(SD,3));
        
        if SD2<=1
            set(handles.tag_std,'backgroundcolor','g');
        elseif SD2>1&SD2<3
            set(handles.tag_std,'backgroundcolor','y');
        else
            set(handles.tag_std,'backgroundcolor','r');
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

% %% FOR AUTOMATIC CORRECTION DURING SR
% rrm = nanmedian(diff(spikes));
% rrmd = diff(spikes,2);
% ii = find(abs(rrmd)>0.66*rrm)
% rr = diff(spikes);
%
% figure
% plot(rr),
% hold on
% plot(ii,rr(ii),'or')

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
    
    sd = diff(handles.spikes);
    SD = nanstd(sd)/handles.ParamSig.frequency*1000;
    sd(sd>nanmedian(sd)+30 | sd<nanmedian(sd)-30)=[];
    SD2 = nanstd(sd)/handles.ParamSig.frequency*1000;
    set(handles.tag_std,'string',num2str(SD,3));
    if SD2<=1
        set(handles.tag_std,'backgroundcolor','g');
    elseif SD2>1&SD2<3
        set(handles.tag_std,'backgroundcolor','y');
    else
        set(handles.tag_std,'backgroundcolor','r');
    end
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


% --- Executes on button press in tag_Eliminate_all.
function tag_Eliminate_all_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Eliminate_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

answ = questdlg('Eliminate all markers currently visualized?');
if isequal(answ,'Yes')
    
    if isfield(handles,'SpikesUpdated');
        SpikesUpdated = handles.SpikesUpdated;
    else
        SpikesUpdated = handles.spikes;
    end
    
    xx = get(handles.tag_ax_sig,'xlim');
    ii_out = SpikesUpdated<xx(2) & SpikesUpdated>xx(1);
    
    
    % % ii_out = nan(1,length(info_struct));
    % for i = 1:length(info_struct)
    %     spikes_out(i) = info_struct(i).Position(1);
    %     ii = find(SpikesUpdated==spikes_out(i));
    %     if ~isempty(ii)
    %         ii_out(i)= min(ii); % min only because in theory ii can be a vector
    %         set(info_struct(i).Target,'visible','off')
    %     end
    % end
    % ii_out(isnan(ii_out))=[];
    plot(handles.tag_ax_sig,SpikesUpdated(ii_out),handles.yym,'xr','linewidth',2,'markersize',10)
    SpikesUpdated(ii_out)=[];
    hold(handles.tag_ax_spikes,'on')
    plot(handles.tag_ax_spikes,SpikesUpdated(2:end),diff(SpikesUpdated),'o-r')
    set(handles.tag_ax_sig,'xtick',SpikesUpdated(2:end),'xticklabel',[],'xgrid','on')
    
    handles.SpikesUpdated=SpikesUpdated;
    
    sd = diff(SpikesUpdated);
    SD = nanstd(sd)/handles.ParamSig.frequency*1000;
    sd(sd>nanmedian(sd)+30 | sd<nanmedian(sd)-30)=[];
    SD2 = nanstd(sd)/handles.ParamSig.frequency*1000;
    set(handles.tag_std,'string',num2str(SD,3));
    if SD2<=1
        set(handles.tag_std,'backgroundcolor','g');
    elseif SD2>1&SD2<3
        set(handles.tag_std,'backgroundcolor','y');
    else
        set(handles.tag_std,'backgroundcolor','r');
    end
    
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



function tag_std_Callback(hObject, eventdata, handles)
% hObject    handle to tag_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_std as text
%        str2double(get(hObject,'String')) returns contents of tag_std as a double


% --- Executes during object creation, after setting all properties.
function tag_std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_control.
function tag_control_Callback(hObject, eventdata, handles)
% hObject    handle to tag_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_control


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
set(handles.tag_std,'string',num2str(SD,3));
if SD2<=1
    set(handles.tag_std,'backgroundcolor','g');
elseif SD2>1&SD2<3
    set(handles.tag_std,'backgroundcolor','y');
else
    set(handles.tag_std,'backgroundcolor','r');
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
ii = handles.spikes<=tt(2)&handles.spikes>=tt(1);
handles.spikes(ii)=[];
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
%
Hwb = waitbar(0,'Estimating Spikes');

%% only interval visualized
signals = handles.sig(ttsamp(1):ttsamp(2),ichan);
if get(handles.tag_do_QRS,'value')==0
    
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
    set(handles.tag_std,'string',num2str(SD,3));
    
    if SD2<=1
        set(handles.tag_std,'backgroundcolor','g');
    elseif SD2>1&SD2<3
        set(handles.tag_std,'backgroundcolor','y');
    else
        set(handles.tag_std,'backgroundcolor','r');
    end
    
    
    % update
    handles.spikes=spikes;
    % end
    handles.SpikesUpdated=spikes;
    guidata(hObject, handles);
    
    %% Sinus Rhythm
else
    
    sp_all = cell(1,size(signals,2));
    sp_all_D = cell(1,size(signals,2));
    x = signals;
    
    waitbar(1/3,Hwb);
    [a,b] = butter(3,[2 30]/(handles.ParamSig.frequency/2));
    xflb = filtfilt(a,b,x);
    
    
    xf = medfilt1(diff(xflb).^2,round(50/1000*handles.ParamSig.frequency));
    %     xf = medfilt1(diff(x).^2,50);
    
    sensthr = str2num(get(handles.tag_sensitivity_thr,'string'));
    h = xf>repmat(max(xf)*sensthr,[size(xf,1) 1]);
    h = imclose(h,ones(round(50/1000*handles.ParamSig.frequency),size(h,2)));
    L = 10;
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
            %             H = ih(j) + [-20:40];
            H = (ih(j)-round(20/1000*handles.ParamSig.frequency)):(eh(j)+round(20/1000*handles.ParamSig.frequency)); % modify to find VT/VF 02/02/2015
            H(H<1 | H>size(xd2,1))=[];
            ii= find([(xd(H(1:length(H)-1),k).*xd(H(2:length(H)),k))<=0 & xd2(H(1:end-1),k)<0])+H(1)+1;
            if length(ii)>1
                %                 [~,id] = min(xd2(ii,k)); % min of second derivative: it
                %                 works for clean and "spiky" QRS
                x2 = detrend(x(H,k));
                [~,id] = max(abs(x2(ii-H(1)))); %
                ii = ii(id);
            end
            if ~isempty(ii)
                sp_all{k}(j) = ii-1;
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
        spikes = spikes/handles.ParamSig.frequency*1000;
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
        set(handles.tag_std,'string',num2str(SD,3));
        
        if SD2<=1
            set(handles.tag_std,'backgroundcolor','g');
        elseif SD2>1&SD2<3
            set(handles.tag_std,'backgroundcolor','y');
        else
            set(handles.tag_std,'backgroundcolor','r');
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
