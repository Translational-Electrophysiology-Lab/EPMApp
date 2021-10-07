function varargout = Alternans_Variability_GUI(varargin)
% ALTERNANS_VARIABILITY_GUI MATLAB code for Alternans_Variability_GUI.fig
%      ALTERNANS_VARIABILITY_GUI, by itself, creates a new ALTERNANS_VARIABILITY_GUI or raises the existing
%      singleton*.
%
%      H = ALTERNANS_VARIABILITY_GUI returns the handle to a new ALTERNANS_VARIABILITY_GUI or the handle to
%      the existing singleton*.
%
%      ALTERNANS_VARIABILITY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALTERNANS_VARIABILITY_GUI.M with the given input arguments.
%
%      ALTERNANS_VARIABILITY_GUI('Property','Value',...) creates a new ALTERNANS_VARIABILITY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Alternans_Variability_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Alternans_Variability_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Alternans_Variability_GUI

% Last Modified by GUIDE v2.5 30-Apr-2014 12:20:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Alternans_Variability_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @Alternans_Variability_GUI_OutputFcn, ...
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


% --- Executes just before Alternans_Variability_GUI is made visible.
function Alternans_Variability_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Alternans_Variability_GUI (see VARARGIN)

% Choose default command line output for Alternans_Variability_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Alternans_Variability_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%%
hObject_main = varargin{1};
CL_pacing = [600:-50:350];
handles_main = guidata(hObject_main); % update
if ~isfield(handles_main,'ALT') % if not already calucated
    % Input from user
    Alternans_Variability_Info(hObject);
    handles = guidata(hObject); % update
    handles.hObject_main = hObject_main;
    set(handles.tag_type,'value',2)
    
    if handles.Alt_input.do_entire_series
        % Initialization
        ALT.ARIw_T=[];ALT.ARIw_SM=[];ALT.STVAR_ARIw=[];
        ALT.ARIa_T=[];ALT.ARIa_SM=[];ALT.STVAR_ARIa=[];
        ALT.RTw_T=[];ALT.RTw_SM=[];ALT.STVAR_RTw=[];
        ALT.RTa_T=[];ALT.RTa_Sm=[];ALT.STVAR_RTa=[];
        ALT.DT_T=[];ALT.DT_SM=[];ALT.STVAR_DT=[];
        ALT.Tp_T=[];ALT.Tp_SM=[];ALT.STVAR_Tp=[];
        %
        if ~handles.Alt_input.do_divide_by_CL
            ALT.ARIw_T_CL=cell(1,length(CL_pacing));ALT.ARIw_SM_CL=cell(1,length(CL_pacing));ALT.STVAR_ARIw_CL=cell(1,length(CL_pacing));
            ALT.ARIa_T_CL=cell(1,length(CL_pacing));ALT.ARIa_SM_CL=cell(1,length(CL_pacing));ALT.STVAR_ARIa_CL=cell(1,length(CL_pacing));
            ALT.RTw_T_CL=cell(1,length(CL_pacing));ALT.RTw_SM_CL=cell(1,length(CL_pacing));ALT.STVAR_RTw_CL=cell(1,length(CL_pacing));
            ALT.RTa_T_CL=cell(1,length(CL_pacing));ALT.RTa_Sm_CL=cell(1,length(CL_pacing));ALT.STVAR_RTa_CL=cell(1,length(CL_pacing));
            ALT.DT_T_CL=cell(1,length(CL_pacing));ALT.DT_SM_CL=cell(1,length(CL_pacing));ALT.STVAR_DT_CL=cell(1,length(CL_pacing));
            ALT.Tp_T_CL=cell(1,length(CL_pacing));ALT.Tp_SM_CL=cell(1,length(CL_pacing));ALT.STVAR_Tp_CL=cell(1,length(CL_pacing));
        end
        totw = handles.Alt_input.do_ARIA+handles.Alt_input.do_RTA+handles.Alt_input.do_DTA+handles.Alt_input.do_ARIA;
        ib = 0;
        wb = waitbar(ib,'Alternans & Variability in Intervals');
        
        % ARI - alternans
        if handles.Alt_input.do_ARIA
            ib=ib+1;waitbar(ib/totw,wb);
            data = handles_main.MarkersC.rt_Wyatt - handles_main.MarkersC.dt;
            [ALT.ARIw_T,ALT.ARIw_SM,ALT.STVAR_ARIw] =  Alternans_Variability_fun(data);
            clear data
            data = handles_main.MarkersC.rt_Alternative - handles_main.MarkersC.dt;
            [ALT.ARIa_T,ALT.ARIa_SM,ALT.STVAR_ARIa] =  Alternans_Variability_fun(data);
            clear data
        end
        % RT - alternans
        if handles.Alt_input.do_RTA
            ib=ib+1;waitbar(ib/totw,wb);
            data = handles_main.MarkersC.rt_Wyatt - repmat( handles_main.spikes(:),[1 size(handles_main.MarkersC.rt_Wyatt,2)]);
            [ALT.RTw_T,ALT.RTw_SM,ALT.STVAR_RTw] =  Alternans_Variability_fun(data);
            clear data
            data = handles_main.MarkersC.rt_Alternative - repmat( handles_main.spikes(:),[1 size(handles_main.MarkersC.rt_Alternative,2)]);
            [ALT.RTa_T,ALT.RTa_SM,ALT.STVAR_RTa] =  Alternans_Variability_fun(data);
            clear data
        end
        
        % AT - alternans
        if handles.Alt_input.do_DTA
            ib=ib+1;waitbar(ib/totw,wb);
            data = handles_main.MarkersC.dt - repmat( handles_main.spikes(:),[1 size(handles_main.MarkersC.dt,2)]);
            [ALT.DT_T,ALT.DT_SM,ALT.STVAR_DT] =  Alternans_Variability_fun(data);
            clear data
        end
        
        % AT - alternans
        if handles.Alt_input.do_TpA
            ib=ib+1;waitbar(ib/totw,wb);
            data = handles_main.MarkersC.tTpeak - repmat( handles_main.spikes(:),[1 size(handles_main.MarkersC.tTpeak,2)]);
            [ALT.Tp_T,ALT.Tp_SM,ALT.STVAR_Tp] =  Alternans_Variability_fun(data);
            clear data
        end
        close(wb)
        
        %% EGM Twave alternans (LLR & SM)
        if isfield(handles_main,'signals_proc');
            Nsurr = [handles.Alt_input.GLRT_n_surr handles.Alt_input.SM_n_surr];
            [ALT.EGMTWA_LLR,ALT.EGMTWA_SM] = EGMTWA_fun(handles_main.signals_proc,handles_main.spikes,handles_main.ParamSig.frequency, handles_main.MarkersC.dt,Nsurr);
        end
        set(handles_main.tag_analysis_special,'backgroundColor', [1 1 1]*.941,'foregroundColor','k')
    end
    
    %% DO Alternans analysis for each Cycle Length
    if handles.Alt_input.do_divide_by_CL
        if ~handles.Alt_input.do_entire_series
            ib = 0;
            wb = waitbar(ib,'Alternans & Variability in Intervals');
        end
        % Initialization ALL
        ALT.ARIw_T_CL=cell(1,length(CL_pacing));ALT.ARIw_SM_CL=cell(1,length(CL_pacing));ALT.STVAR_ARIw_CL=cell(1,length(CL_pacing));
        ALT.ARIa_T_CL=cell(1,length(CL_pacing));ALT.ARIa_SM_CL=cell(1,length(CL_pacing));ALT.STVAR_ARIa_CL=cell(1,length(CL_pacing));
        ALT.RTw_T_CL=cell(1,length(CL_pacing));ALT.RTw_SM_CL=cell(1,length(CL_pacing));ALT.STVAR_RTw_CL=cell(1,length(CL_pacing));
        ALT.RTa_T_CL=cell(1,length(CL_pacing));ALT.RTa_Sm_CL=cell(1,length(CL_pacing));ALT.STVAR_RTa_CL=cell(1,length(CL_pacing));
        ALT.DT_T_CL=cell(1,length(CL_pacing));ALT.DT_SM_CL=cell(1,length(CL_pacing));ALT.STVAR_DT_CL=cell(1,length(CL_pacing));
        ALT.Tp_T_CL=cell(1,length(CL_pacing));ALT.Tp_SM_CL=cell(1,length(CL_pacing));ALT.STVAR_Tp_CL=cell(1,length(CL_pacing));
        for iCL = 1:length(CL_pacing)
            wb = waitbar(ib,['Alternans & Variability in Intervals [CL=',num2str(CL_pacing(iCL)),']']);
            
            % select data with desired CL only
            iihb = find(diff(handles_main.spikes)>CL_pacing(iCL)-10 & diff(handles_main.spikes)<CL_pacing(iCL)+10);
            if length(iihb)>10 % at least 10 heart beats
                if diff(iihb)>1
                    error('[CL=',num2str(CL_pacing(iCL)),': not continuous sequence]')
                end
                
                % Initialization
                ALT.ARIw_T_CL{iCL}=[];ALT.ARIw_SM_CL{iCL}=[];ALT.STVAR_ARIw_CL{iCL}=[];
                ALT.ARIa_T_CL{iCL}=[];ALT.ARIa_SM_CL{iCL}=[];ALT.STVAR_ARIa_CL{iCL}=[];
                ALT.RTw_T_CL{iCL}=[];ALT.RTw_SM_CL{iCL}=[];ALT.STVAR_RTw_CL{iCL}=[];
                ALT.RTa_T_CL{iCL}=[];ALT.RTa_Sm_CL{iCL}=[];ALT.STVAR_RTa_CL{iCL}=[];
                ALT.DT_T_CL{iCL}=[];ALT.DT_SM_CL{iCL}=[];ALT.STVAR_DT_CL{iCL}=[];
                ALT.Tp_T_CL{iCL}=[];ALT.Tp_SM_CL{iCL}=[];ALT.STVAR_Tp_CL{iCL}=[];
                
                if ~handles.Alt_input.do_entire_series
                    ALT.ARIw_T=[];ALT.ARIw_SM=[];ALT.STVAR_ARIw=[];
                    ALT.ARIa_T=[];ALT.ARIa_SM=[];ALT.STVAR_ARIa=[];
                    ALT.RTw_T=[];ALT.RTw_SM=[];ALT.STVAR_RTw=[];
                    ALT.RTa_T=[];ALT.RTa_Sm=[];ALT.STVAR_RTa=[];
                    ALT.DT_T=[];ALT.DT_SM=[];ALT.STVAR_DT=[];
                    ALT.Tp_T=[];ALT.Tp_SM=[];ALT.STVAR_Tp=[];
                end
                %
                totw = handles.Alt_input.do_ARIA+handles.Alt_input.do_RTA+handles.Alt_input.do_DTA+handles.Alt_input.do_ARIA;
                ib = 0;
                
                % ARI - alternans
                if handles.Alt_input.do_ARIA
                    ib=ib+1;waitbar(ib/totw,wb);
                    data = handles_main.MarkersC.rt_Wyatt - handles_main.MarkersC.dt;
                    data = data(iihb,:); % only data with desired CL
                    [ALT.ARIw_T_CL{iCL},ALT.ARIw_SM_CL{iCL},ALT.STVAR_ARIw_CL{iCL}] =  Alternans_Variability_fun(data);
                    clear data
                    data = handles_main.MarkersC.rt_Alternative - handles_main.MarkersC.dt;
                    data = data(iihb,:); % only data with desired CL
                    [ALT.ARIa_T_CL{iCL},ALT.ARIa_SM_CL{iCL},ALT.STVAR_ARIa_CL{iCL}] =  Alternans_Variability_fun(data);
                    clear data
                end
                % RT - alternans
                if handles.Alt_input.do_RTA
                    ib=ib+1;waitbar(ib/totw,wb);
                    data = handles_main.MarkersC.rt_Wyatt - repmat( handles_main.spikes(:),[1 size(handles_main.MarkersC.rt_Wyatt,2)]);
                    data = data(iihb,:); % only data with desired CL
                    [ALT.RTw_T_CL{iCL},ALT.RTw_SM_CL{iCL},ALT.STVAR_RTw_CL{iCL}] =  Alternans_Variability_fun(data);
                    clear data
                    data = handles_main.MarkersC.rt_Alternative - repmat( handles_main.spikes(:),[1 size(handles_main.MarkersC.rt_Alternative,2)]);
                    data = data(iihb,:); % only data with desired CL
                    [ALT.RTa_T_CL{iCL},ALT.RTa_SM_CL{iCL},ALT.STVAR_RTa_CL{iCL}] =  Alternans_Variability_fun(data);
                    clear data
                end
                
                % AT - alternans
                if handles.Alt_input.do_DTA
                    ib=ib+1;waitbar(ib/totw,wb);
                    data = handles_main.MarkersC.dt - repmat( handles_main.spikes(:),[1 size(handles_main.MarkersC.dt,2)]);
                    data = data(iihb,:); % only data with desired CL
                    [ALT.DT_T_CL{iCL},ALT.DT_SM_CL{iCL},ALT.STVAR_DT_CL{iCL}] =  Alternans_Variability_fun(data);
                    clear data
                end
                
                % T-peak - alternans
                if handles.Alt_input.do_TpA
                    ib=ib+1;waitbar(ib/totw,wb);
                    data = handles_main.MarkersC.tTpeak - repmat( handles_main.spikes(:),[1 size(handles_main.MarkersC.tTpeak,2)]);
                    data = data(iihb,:); % only data with desired CL
                    [ALT.Tp_T_CL{iCL},ALT.Tp_SM_CL{iCL},ALT.STVAR_Tp_CL{iCL}] =  Alternans_Variability_fun(data);
                    clear data
                end
            end
            close(wb)
            
            %% EGM Twave alternans (LLR & SM)
            if isfield(handles_main,'signals_proc')
                Nsurr = [handles.Alt_input.GLRT_n_surr handles.Alt_input.SM_n_surr];
                [ALT.EGMTWA_LLR,ALT.EGMTWA_SM] = EGMTWA_fun(handles_main.signals_proc,handles_main.spikes,handles_main.ParamSig.frequency, handles_main.MarkersC.dt,Nsurr);
            end
            set(handles_main.tag_analysis_special,'backgroundColor', [1 1 1]*.941,'foregroundColor','k')
        end
    end
    
    
    
    %% Labels
    for i= 1:length(handles_main.ParamSig.Label)
        lab(i) = {['[',num2str(i),'] ',handles_main.ParamSig.Label{i},' / ',num2str(handles_main.SNR(i),3),'dB']};
    end
    set(handles.tag_List,'string',lab)
    set(handles.tag_ax_marker,'xtick',[],'ytick',[])
    set(handles.tag_ax_sig,'xtick',[],'ytick',[])
    %%
    guidata(hObject,handles)
    handles_main.ALT =  ALT;
    guidata(hObject_main,handles_main)
    %%
    if isempty(handles_main.ALT.RTw_T)
        set(handles.tag_plot_RTw,'foregroundColor',[.7 .7 .7])
    end
    if isempty(handles_main.ALT.ARIw_T)
        set(handles.tag_plot_ARIw,'foregroundColor',[.7 .7 .7])
    end
    if isempty(handles_main.ALT.ARIa_T)
        set(handles.tag_plot_ARIa,'foregroundColor',[.7 .7 .7])
    end
    
    if isempty(handles_main.ALT.DT_T)
        set(handles.tag_plot_DT,'foregroundColor',[.7 .7 .7])
    end
    if isempty(handles_main.ALT.Tp_T)
        set(handles.tag_plot_Tp,'foregroundColor',[.7 .7 .7])
    end
    
else
    handles.hObject_main = hObject_main;
    guidata(hObject,handles)
    for i= 1:length(handles_main.ParamSig.Label)
        lab(i) = {[handles_main.ParamSig.Label{i},' / ',num2str(handles_main.SNR(i),3),'dB']};
    end
    set(handles.tag_List,'string',lab)
    set(handles.tag_ax_marker,'xtick',[],'ytick',[])
    set(handles.tag_ax_sig,'xtick',[],'ytick',[])
end

% --- Outputs from this function are returned to the command line.
function varargout = Alternans_Variability_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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


% --- Executes on selection change in tag_List.
function tag_List_Callback(hObject, eventdata, handles)
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


% --- Executes on button press in tag_plot.
function tag_plot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.tag_type,'value')==1
    set(handles.tag_type,'value',2)
end
data = guidata(handles.hObject_main);
ic = get(handles.tag_List,'value');

if get(handles.tag_filt,'value')
    sig = data.signals_proc(:,ic)./repmat(sqrt(nansum(data.signals_proc(:,ic).^2)),[size(data.signals_proc,1) 1]);
else
    sig = data.signals(:,ic)./repmat(sqrt(nansum(data.signals(:,ic).^2)),[size(data.signals,1) 1]);
end

hold(handles.tag_ax_sig,'off');
aasig = plot(handles.tag_ax_sig,[1 : size(sig,1)]/data.ParamSig.frequency*1000,sig); % [ms]
set(handles.tag_ax_sig,'ytick',[])
legend(aasig,data.ParamSig.Label{ic})
hold(handles.tag_ax_marker,'off')

nB = str2double(get(handles.tag_n_beats_Alt,'string'))-5;
if get(handles.tag_type,'value')==2
    %% ari wyatt
    if get(handles.tag_plot_ARIw,'value')&~isempty(data.ALT.ARIw_T)
        aaari=plot(handles.tag_ax_marker,data.MarkersC.rt_Wyatt(:,ic),data.MarkersC.rt_Wyatt(:,ic)-data.MarkersC.dt(:,ic),'*--');
        
        iirt = round(data.MarkersC.rt_Wyatt(:,ic)/1000*data.ParamSig.frequency);
        hold(handles.tag_ax_sig,'on')
        for i=1:length(ic)
            plot(handles.tag_ax_sig, data.MarkersC.rt_Wyatt(~isnan(iirt(:,i)),ic(i)), sig(iirt(~isnan(iirt(:,i)),i),i),'x','color',get(aaari(i),'color'),'linewidth',2,'markersize',10);
        end
        clear iirt
        iidt = round(data.MarkersC.dt(:,ic)/1000*data.ParamSig.frequency);
        for i=1:length(ic)
            plot(handles.tag_ax_sig, data.MarkersC.dt(~isnan(iidt(:,i)),ic(i)), sig(iidt(~isnan(iidt(:,i)),i),i),'o','color',get(aaari(i),'color'),'markerfacecolor',get(aaari(i),'color'),'linewidth',2);
        end
        clear iidt
        
        ii = data.ALT.ARIw_T.Info_beats_pos(ic,nB);
        hold(handles.tag_ax_marker,'on')
        for i = 1:length(ii)
            if ~isempty(ii{i})
                plot(handles.tag_ax_marker,data.MarkersC.rt_Wyatt(ii{i},ic(i)),data.MarkersC.rt_Wyatt(ii{i},ic(i))-data.MarkersC.dt(ii{i},ic(i)),'*','linewidth',2,'markersize',10,'color',get(aaari(i),'color'))
                plot(handles.tag_ax_marker,data.MarkersC.rt_Wyatt(ii{i},ic(i)),data.MarkersC.rt_Wyatt(ii{i},ic(i))-data.MarkersC.dt(ii{i},ic(i)),'o','linewidth',.5,'markersize',10,'color','k')
            end
        end
    end
    %% ari alternative
    if get(handles.tag_plot_ARIa,'value')&~isempty(data.ALT.ARIa_T)
        aaari=plot(handles.tag_ax_marker,data.MarkersC.rt_Alternative(:,ic),data.MarkersC.rt_Alternative(:,ic)-data.MarkersC.dt(:,ic),'*--');
        
        iirt = round(data.MarkersC.rt_Alternative(:,ic)/1000*data.ParamSig.frequency);
        hold(handles.tag_ax_sig,'on')
        for i=1:length(ic)
            plot(handles.tag_ax_sig, data.MarkersC.rt_Alternative(~isnan(iirt(:,i)),ic(i)), sig(iirt(~isnan(iirt(:,i)),i),i),'x','color',get(aaari(i),'color'),'linewidth',2,'markersize',10);
        end
        clear iirt
        iidt = round(data.MarkersC.dt(:,ic)/1000*data.ParamSig.frequency);
        for i=1:length(ic)
            plot(handles.tag_ax_sig, data.MarkersC.dt(~isnan(iidt(:,i)),ic(i)), sig(iidt(~isnan(iidt(:,i)),i),i),'o','color',get(aaari(i),'color'),'markerfacecolor',get(aaari(i),'color'),'linewidth',2);
        end
        clear iidt
        
        ii = data.ALT.ARIa_T.Info_beats_pos(ic,nB);
        hold(handles.tag_ax_marker,'on')
        for i = 1:length(ii)
            if ~isempty(ii{i})
                plot(handles.tag_ax_marker,data.MarkersC.rt_Alternative(ii{i},ic(i)),data.MarkersC.rt_Alternative(ii{i},ic(i))-data.MarkersC.dt(ii{i},ic(i)),'*','linewidth',2,'markersize',10,'color',get(aaari(i),'color'))
                plot(handles.tag_ax_marker,data.MarkersC.rt_Alternative(ii{i},ic(i)),data.MarkersC.rt_Alternative(ii{i},ic(i))-data.MarkersC.dt(ii{i},ic(i)),'o','linewidth',.5,'markersize',10,'color','k')
            end
        end
    end
    %%
    if get(handles.tag_plot_RTw,'value')&~isempty(data.ALT.RTw_T)
        aart=plot(handles.tag_ax_marker,data.MarkersC.rt_Wyatt(:,ic),data.MarkersC.rt_Wyatt(:,ic)- repmat(data.spikes(:),[1,length(ic)]),'x--');
        ii = data.ALT.RTw_T.Info_beats_pos(ic,nB);
        hold(handles.tag_ax_marker,'on')
        for i = 1:length(ii)
            if ~isempty(ii{i})
                plot(handles.tag_ax_marker,data.MarkersC.rt_Wyatt(ii{i},ic(i)),data.MarkersC.rt_Wyatt(ii{i},ic(i))- data.spikes(ii{i})','x','linewidth',2,'markersize',10,'color',get(aart(i),'color'))
                plot(handles.tag_ax_marker,data.MarkersC.rt_Wyatt(ii{i},ic(i)),data.MarkersC.rt_Wyatt(ii{i},ic(i))- data.spikes(ii{i})','o','linewidth',.5,'markersize',10,'color','k')
                
            end
        end
    end
    %%
    if get(handles.tag_plot_DT,'value')&~isempty(data.ALT.DT_T)
        aadt = plot(handles.tag_ax_marker,data.MarkersC.dt(:,ic),data.MarkersC.dt(:,ic)- repmat(data.spikes(:),[1,length(ic)]),'.--');
        ii = data.ALT.DT_T.Info_beats_pos(ic,nB);
        hold(handles.tag_ax_marker,'on')
        for i = 1:length(ii)
            if ~isempty(ii{i})
                plot(handles.tag_ax_marker,data.MarkersC.dt(ii{i},ic(i)),data.MarkersC.dt(ii{i},ic(i))- data.spikes(ii{i})','o','markersize',10)
                plot(handles.tag_ax_marker,data.MarkersC.dt(ii{i},ic(i)),data.MarkersC.dt(ii{i},ic(i))-data.spikes(ii{i})','o','linewidth',.5,'markersize',10,'color','k')
            end
        end
    end
    %%
    if get(handles.tag_plot_Tp,'value')&~isempty(data.ALT.Tp_T)
        aatp = plot(handles.tag_ax_marker,data.MarkersC.tTpeak(:,ic),data.MarkersC.tTpeak(:,ic)- repmat(data.spikes(:),[1,length(ic)]),'v--');
        ii = data.ALT.Tp_T.Info_beats_pos(ic,nB);
        hold(handles.tag_ax_marker,'on')
        for i = 1:length(ii)
            if ~isempty(ii{i})
                plot(handles.tag_ax_marker,data.MarkersC.tTpeak(ii{i},ic(i)),data.MarkersC.tTpeak(ii{i},ic(i))- data.spikes(ii{i})','v','markeredgecolor',get(aatp(i),'color'),'markersize',10)
                plot(handles.tag_ax_marker,data.MarkersC.tTpeak(ii{i},ic(i)),data.MarkersC.tTpeak(ii{i},ic(i))-data.spikes(ii{i})','o','linewidth',.5,'markersize',10,'color','k')
            end
        end
    end
    %%
    %% write summary
    summ = cell(1);
    for i = 1:length(ic)
        ARIaltdata_W = eval(['data.ALT.ARIw_T.nbeats_',num2str(nB+5)]);
        ARIaltdata_A = eval(['data.ALT.ARIa_T.nbeats_',num2str(nB+5)]);
        DTaltdata = eval(['data.ALT.DT_T.nbeats_',num2str(nB+5)]);
        
        str = {['--- IC ',data.ParamSig.Label{ic(i)},' ---'],...
            [' - ARI(W)-Alt mean=',num2str(ARIaltdata_W.mean(ic(i)),2),' ms'],...
            [' - ARI(W)-Alt med=',num2str(ARIaltdata_W.median(ic(i)),2),' ms'],...
            [' - ARI(W)-Alt #beats=',num2str(ARIaltdata_W.nbeats(ic(i)),2)],...
            [' - ARI(W)-SM K-score=',num2str(data.ALT.ARIw_SM.K(ic(i)),2),'/',num2str(data.ALT.ARIw_SM.Kth(2),2)],...
            [' - ARI(W)-Var mean=',num2str( data.ALT.STVAR_ARIw.mean(ic(i)),2),' ms'],...
            [' - ARI(W)-Var med=',num2str(data.ALT.STVAR_ARIw.median(ic(i)),2),' ms'],' ',...
            [' - ARI(A)-Alt mean=',num2str(ARIaltdata_A.mean(ic(i)),2),' ms'],...
            [' - ARI(A)-Alt med=',num2str(ARIaltdata_A.median(ic(i)),2),' ms'],...
            [' - ARI(A)-Alt #beats=',num2str(ARIaltdata_A.nbeats(ic(i)),2)],...
            [' - ARI(A)-SM K-score=',num2str(data.ALT.ARIa_SM.K(ic(i)),2),'/',num2str(data.ALT.ARIa_SM.Kth(2),2)],...
            [' - ARI(A)-Var mean=',num2str( data.ALT.STVAR_ARIa.mean(ic(i)),2),' ms'],...
            [' - ARI(A)-Var med=',num2str(data.ALT.STVAR_ARIa.median(ic(i)),2),' ms'],' ',...
            [' - DT-Alt mean=',num2str(DTaltdata.mean(ic(i)),2),' ms'],...
            [' - DT-Alt med=',num2str(DTaltdata.median(ic(i)),2),' ms'],...
            [' - DT-Alt #beats=',num2str(DTaltdata.nbeats(ic(i)),2)],...
            [' - DT-SM K-score=',num2str(data.ALT.DT_SM.K(ic(i)),2),'/',num2str(data.ALT.DT_SM.Kth(2),2)],...
            [' - DT-Var mean=',num2str( data.ALT.STVAR_DT.mean(ic(i)),2),' ms'],...
            [' - DT-Var med=',num2str(data.ALT.STVAR_DT.median(ic(i)),2),' ms'],' '};
        if ~isempty(data.ALT.EGMTWA_SM.ks_th)
            str2 = [' - EGM-TWA-SM K-score=',num2str(data.ALT.EGMTWA_SM.ks(ic(i)),2),'/',num2str(data.ALT.EGMTWA_SM.ks_th(2),2)];
        else
            str2 = [' - EGM-TWA-SM K-score=',num2str(data.ALT.EGMTWA_SM.ks(ic(i)),2),'/NP'];
        end
        
        if ~isempty(data.ALT.EGMTWA_LLR.z_th)
            str3 = [' - EGM-TWA-LLR z=',num2str(data.ALT.EGMTWA_LLR.z(ic(i)),2),'/',num2str(data.ALT.EGMTWA_LLR.z_th(2),2)];
        else
            str3 = [' - EGM-TWA-SM K-score=',num2str(data.ALT.EGMTWA_LLR.z(ic(i)),2),'/NP'];
        end
        
        summ = {summ{:},str{:},str2,str3,' '};
    end
    set(handles.tag_summary,'string',summ)
end
if get(handles.tag_type,'value')==3
    hold(handles.tag_ax_marker,'off')
    if get(handles.tag_plot_ARIw,'value')
        arivar = abs(diff(data.MarkersC.rt_Wyatt(:,ic)-data.MarkersC.dt(:,ic)));
        aaari=plot(handles.tag_ax_marker,data.MarkersC.rt_Wyatt(2:end,ic),arivar,'*--');
        hold(handles.tag_ax_marker,'on')
        plot(handles.tag_ax_marker,get(handles.tag_ax_sig,'xlim'),nanmean(arivar)'*[1 1],'linewidth',3)
    end
end

linkaxes([handles.tag_ax_sig,handles.tag_ax_marker],'x')
% --- Executes on selection change in tag_type.
function tag_type_Callback(hObject, eventdata, handles)
% hObject    handle to tag_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_type


% --- Executes during object creation, after setting all properties.
function tag_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_plot_ARIw.
function tag_plot_ARIw_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_ARIw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_ARIw


% --- Executes on button press in tag_plot_DT.
function tag_plot_DT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_DT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_DT


% --- Executes on button press in tag_plot_RTw.
function tag_plot_RTw_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_RTw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_RTw


% --- Executes on button press in tag_plot_Tp.
function tag_plot_Tp_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_Tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_Tp



function tag_summary_Callback(hObject, eventdata, handles)
% hObject    handle to tag_summary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_summary as text
%        str2double(get(hObject,'String')) returns contents of tag_summary as a double


% --- Executes during object creation, after setting all properties.
function tag_summary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_summary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_n_beats_Alt_Callback(hObject, eventdata, handles)
% hObject    handle to tag_n_beats_Alt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_n_beats_Alt as text
%        str2double(get(hObject,'String')) returns contents of tag_n_beats_Alt as a double


% --- Executes during object creation, after setting all properties.
function tag_n_beats_Alt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_n_beats_Alt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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


% --- Executes on button press in tag_all_channels.
function tag_all_channels_Callback(hObject, eventdata, handles)
% hObject    handle to tag_all_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = guidata(handles.hObject_main);

nb = str2double(get(handles.tag_n_beats_Alt,'string'));
ppT = cell(1,6);
ppSM = cell(1,6);
ppSM_thr = cell(1,6);
titles_ax = {'ARI (Wyatt)','ARI (Alternative)','RT (Wyatt)','RT (Alternative)','DT','Tpe'};
if get(handles.tag_plot_ARIw,'value')
    ppT{1} = eval(['data.ALT.ARIw_T.nbeats_',num2str(nb)]);
    ppSM{1} = data.ALT.ARIw_SM.K;
    ppSM_thr{1} = data.ALT.ARIw_SM.Kth;
end
if get(handles.tag_plot_ARIa,'value')
    ppT{2} = eval(['data.ALT.ARIa_T.nbeats_',num2str(nb)]);
    ppSM{2} = data.ALT.ARIa_SM.K;
    ppSM_thr{2} = data.ALT.ARIa_SM.Kth;
end
if get(handles.tag_plot_RTw,'value')
    ppT{3} = eval(['data.ALT.RTw_T.nbeats_',num2str(nb)]);
    ppSM{3} = data.ALT.RTw_SM.K;
    ppSM_thr{3} = data.ALT.RTw_SM.Kth;
end
if get(handles.tag_plot_RTa,'value')
    ppT{4} = eval(['data.ALT.RTw_T.nbeats_',num2str(nb)]);
    ppSM{4} = data.ALT.RTa_SM.K;
    ppSM_thr{4} = data.ALT.RTa_SM.Kth;
end
if get(handles.tag_plot_DT,'value')
    ppT{5} = eval(['data.ALT.DT_T.nbeats_',num2str(nb)]);
    ppSM{5} = data.ALT.DT_SM.K;
    ppSM_thr{5} = data.ALT.DT_SM.Kth;
end
if get(handles.tag_plot_Tp,'value')
    ppT{6} = eval(['data.ALT.Tp_T.nbeats_',num2str(nb)]);
    ppSM{6} = data.ALT.Tp_SM.K;
    ppSM_thr{6} = data.ALT.Tp_SM.Kth;
end

iiok = find(data.SNR>15);

% plot Alternans T (median alternans in sequences of nb beats)
for i = 1:length(ppT)
    if ~isempty(ppT{i})
        figure
        ax(1) = subplot(211);
        a1 = plot([1:length(ppT{i}.median)],ppT{i}.median);
        hold on
        a12 = plot(iiok,ppT{i}.median(iiok));
        % --
        ax(2) = subplot(212);
        a2 = plot([1:length(ppT{i}.nbeats)],ppT{i}.nbeats);
        hold on
        a22 = plot(iiok,ppT{i}.nbeats(iiok));
        
        set([a1,a2],'marker','o','color',[.8 .8 .8],'markerfacecolor',[.8 .8 .8],'markersize',5,'linestyle','--');
        set([a12,a22],'marker','o','color','b','markerfacecolor','b','markersize',6,'linestyle','none');
        legend([a12,a1],'SNR>15','SNR<15')
        legend([a22,a2],'SNR>15','SNR<15')
        
        ylabel(ax(1),'Alternans (median)')
        ylabel(ax(2),'beats # with alternans')
        xlabel(ax(2),'channles')
        linkaxes(ax,'x')
        title(ax(1),[titles_ax{i},' Alternans [Temporal, in ',num2str(nb),' beats]'])
        clear ax a1* a2*
    end
end

% plot Alternans SM
for i = 1:length(ppT)
    if ~isempty(ppSM{i})
        figure
        a1 = plot([1:length(ppSM{i})],ppSM{i});
        hold on
        a12 = plot(iiok,ppSM{i}(iiok));
        if ~isempty(ppSM_thr{i})
            a13 = plot([1 length(ppSM{i})],[1 1]*ppSM_thr{i}(2),'r')
            legend([a12,a1,a13],'SNR>15','SNR<15','thr-5%')
        else
            legend([a12,a1],'SNR>15','SNR<15')
        end
        set([a1],'marker','o','color',[.8 .8 .8],'markerfacecolor',[.8 .8 .8],'markersize',5,'linestyle','--');
        set([a12],'marker','o','color','b','markerfacecolor','b','markersize',6,'linestyle','none');
        
        
        ylabel(gca,'Alternans (K-score SM)')
        xlabel(gca,'channles')
        title(gca,[titles_ax{i},' Alternans [Spectral Method]'])
    end
end

%% PLOT EGM-TWA LLR
clear ax a1* a2*
figure
ax(1)=subplot(211);
a1 = plot(data.ALT.EGMTWA_LLR.v);
hold on
a12 = plot(iiok,data.ALT.EGMTWA_LLR.v(iiok));
if ~isempty(data.ALT.EGMTWA_LLR.v_th)
    a13 = plot(data.ALT.EGMTWA_LLR.v_th(:,2),'r');
else
    a13=[];
end
ax(2)=subplot(212);
a2 = plot(data.ALT.EGMTWA_LLR.z);
hold on
a22 = plot(iiok,data.ALT.EGMTWA_LLR.z(iiok));
if ~isempty(data.ALT.EGMTWA_LLR.v_th)
    a23 = plot(data.ALT.EGMTWA_LLR.z_th(:,2),'r');
else
    a23=[];
end


set([a1,a2],'marker','o','color',[.8 .8 .8],'markerfacecolor',[.8 .8 .8],'markersize',5,'linestyle','--');
set([a12,a22],'marker','o','color','b','markerfacecolor','b','markersize',6,'linestyle','none');

legend([a12,a1,a13],'SNR>15','SNR<15','thr-5%')
legend([a22,a2,a23],'SNR>15','SNR<15','thr-5%')


ylabel(ax(1),'Alternans (LLR v)')
ylabel(ax(2),'Alternans (LLR z)')
xlabel(ax(2),'channles')
linkaxes(ax,'x')
title(ax(1),[' EGM-Twave Alternans [LLR]'])
clear ax a1* a2*

%% PLOT EGM-TWA SM
clear ax a1* a2*
figure
ax(1)=subplot(211);
a1 = plot(data.ALT.EGMTWA_SM.Ptwa);
hold on
a12 = plot(iiok,data.ALT.EGMTWA_SM.Ptwa(iiok));
if ~isempty( data.ALT.EGMTWA_SM.Ptwa_th)
    a13 = plot(data.ALT.EGMTWA_SM.Ptwa_th(:,2),'r');
else
    a13 = [];
end
ax(2)=subplot(212);
a2 = plot(data.ALT.EGMTWA_SM.ks);
hold on
a22 = plot(iiok,data.ALT.EGMTWA_SM.ks(iiok));
if ~isempty( data.ALT.EGMTWA_SM.ks_th)
    a23 = plot(data.ALT.EGMTWA_SM.ks_th(:,2),'r');
else
    a23 = [];
end





set([a1,a2],'marker','o','color',[.8 .8 .8],'markerfacecolor',[.8 .8 .8],'markersize',5,'linestyle','--');
set([a12,a22],'marker','o','color','b','markerfacecolor','b','markersize',6,'linestyle','none');
legend([a12,a1,a13],'SNR>15','SNR<15','thr-5%')
legend([a22,a2,a23],'SNR>15','SNR<15','thr-5%')


ylabel(ax(1),'Alternans (SM v)')
ylabel(ax(2),'Alternans (SM z)')
xlabel(ax(2),'channles')
linkaxes(ax,'x')
title(ax(1),['EGM-Twave Alternans [SM]'])
clear ax a1* a2*
%%



% --- Executes on button press in tag_done.
function tag_done_Callback(hObject, eventdata, handles)
% hObject    handle to tag_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1)

% --- Executes on button press in tag_plot_ARIa.
function tag_plot_ARIa_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_ARIa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_ARIa


% --- Executes on button press in tag_plot_RTa.
function tag_plot_RTa_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_RTa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_RTa
