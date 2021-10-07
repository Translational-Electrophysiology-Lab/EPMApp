function varargout = EPMApp(varargin)
% EPMAPP MATLAB code for EPMApp.fig
%      EPMAPP, by itself, creates a new EPMAPP or raises the existing
%      singleton*.
%
%      H = EPMAPP returns the handle to a new EPMAPP or the handle to
%      the existing singleton*.
%
%      EPMAPP('CALLBACK',hObjecguidet,eventData,handles,...) calls the local
%      function named CALLBACK in EPMAPP.M with the given input arguments.
%
%      EPMAPP('Property','Value',...) creates a new EPMAPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EPMApp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_negm_main_v4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".0
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EPMApp

% Last Modified by GUIDE v2.5 21-Sep-2021 12:15:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @EPMApp_OpeningFcn, ...
    'gui_OutputFcn',  @EPMApp_OutputFcn, ...
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


% --- Executes just before EPMApp is made visible.
function EPMApp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EPMApp (see VARARGIN)

% Choose default command line output for EPMApp
handles.output = hObject;
handles.LOG =  {datestr(now)};
handles.outputname =  [];
addpath([pwd,filesep,'GUI_egm_mFiles'])
addpath([pwd,filesep,'GUIs'])
guidata(hObject, handles);
%
% % UIWAIT makes EPMApp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EPMApp_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;



function tag_no_chan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_no_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of tag_no_chan as text
%        str2double(get(hObject,'String')) returns contents of tag_no_chan as a double


% --- Executes during object creation, after setting all properties.
function tag_no_chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_no_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in tag_new_marker.
function tag_new_marker_Callback(hObject, eventdata, handles)
% hObject    handle to tag_new_marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'UpdateFcn',[])
% pause




% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject)
else
    delete(hObject);
end



% --- Executes on button press in tag_import_data.
function tag_import_data_Callback(hObject, eventdata, handles)
% hObject    handle to tag_import_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% clear data from previous analysis
v = fieldnames(handles);
ii = find(~strncmp(v,'tag',3)&~strncmp(v,'fig',3)&~strncmp(v,'text',4)&~strncmp(v,'uip',3)&~strcmp(v,'output')&~strcmp(v,'LOG'));
for i=1:length(ii)
    handles = rmfield(handles,v{ii(i)});
end
clear i ii
%%
% cc1 = get(handles.tag_import_data,'BackgroundColor');
% cc2 = get(handles.tag_import_data,'fontweight');
set(handles.tag_import_data,'BackgroundColor',[1 0 0],'fontweight','bold')

format = handles.tag_import_type.String{handles.tag_import_type.Value};
switch format
    case 'Bard(*.txt)'
        [fnamesig,pname] = uigetfile('*.txt', 'Choose a Signal File');
        filename = [pname fnamesig];
        %
        if sum(filename)==0
            set(handles.tag_import_data,'BackgroundColor',[1 1 1]*.941,'fontweight','normal')
            handles.LOG = [handles.LOG,'- No file loaded'];
            set(handles.tag_Log,'string',handles.LOG)
            return
        end
        %
        handles.LOG = [handles.LOG,['- Loading ',filename,],[' ... please wait ...']];
        set(handles.tag_Log,'string',handles.LOG)
        %
        [signals_raw,ParamSig] = LoadBardSignals_moCathLab(filename);
        
        handles.LOG = [handles.LOG,['- Loaded: ',fnamesig]];
        set(handles.tag_Log,'string',handles.LOG)
        
        
        st = {'cath-lab','mo_sock2','mo_sock1','new_sock4','old_sock4','old_sock6','old_sock6_newPCB','new_sock','new_sock4_Stefan','mo_sock3'};
        v = listdlg('ListString',st,'SelectionMode','single','PromptString','Select a geometry');
        handles.geoname = st{v};
        ParamSig.geoname = st{v};
        %
        handles.LOG = [handles.LOG,['- selected geometry:',st{v}]];
        set(handles.tag_Log,'string',handles.LOG)
        %
        %     set(handles.tag_name,'string',['- ',fnamesig]);
        
        ParamSig.name_file = fnamesig;
        ParamSig.name_dir = pname;
        
        
        handles.signals_raw=signals_raw;
        clear signals_raw
        
    case 'Matlab(*.mat)'
        
        hh = helpdlg({'Imported file should contain:','signals_raw : matrix TxNc, with Nc # of channel','ParamSig : struct as returned by LoadBardSignals_moCathLab.m'},'HELP');
        
        [fnamesig,pname] = uigetfile('*.mat', 'Choose a Signal File');
        filename = [pname fnamesig];
        
        %
        if sum(filename)==0
            set(handles.tag_import_data,'BackgroundColor',[1 1 1]*.941,'fontweight','normal')
            handles.LOG = [handles.LOG,'- No file loaded'];
            set(handles.tag_Log,'string',handles.LOG)
            close(hh)
            return
        end
        %
        
        load(filename)
        %
        handles.LOG = [handles.LOG,['- Loading ',filename,':']];
        set(handles.tag_Log,'string',handles.LOG)
        
        %
        if exist('geoname','var')
            st=geoname;
        else
            st = 'Unknown';
        end
        set(handles.tag_name,'string',['- ',fnamesig]);
        %
        V = whos('-file',[filename]);
        for i=1:length(V)
            Vfields{i} = V(i).name;
            handles.(Vfields{i}) = eval(Vfields{i});
        end
        %
        handles.LOG = [handles.LOG,[' -------- ',sort(Vfields),' -------- ']];
        set(handles.tag_Log,'string',handles.LOG)
        %
        
        if exist('ParamSig','var')
            if ~isfield(ParamSig,'Label')
                %             ParamSig.Label = [{'ALL'},num2cell(1:size(signals_raw,2))];
                for kk = 1:size(signals,2)
                    ParamSig.Label{kk} = ['chan-',num2str(kk)];
                end
            end
            if ~isfield(ParamSig,'name_file')
                ParamSig.name_file = fnamesig;
                ParamSig.name_dir = pname;
            end
            if ~isfield(ParamSig,'frequency')
                answ = inputdlg('Enter sampling frequency [Hz]','sampling frequency',1,{'1000'});
                ParamSig.frequency = str2double(answ{1});
            end
        else
            ParamSig.name_file = fnamesig;
            ParamSig.name_dir = pname;
            answ = inputdlg('Enter sampling frequency [Hz]','sampling frequency',1,{'1000'});
            ParamSig.frequency = str2double(answ{1});
        end
        clear kk
        close(hh)
        
        ParamSig.geoname = st;
        %
    case 'Ensite(*.csv)'
        [fnamesig,pname] = uigetfile('*.csv', 'Choose a Signal File');
        filename = [pname fnamesig];
        %
        if sum(filename)==0
            set(handles.tag_import_data,'BackgroundColor',[1 1 1]*.941,'fontweight','normal')
            handles.LOG = [handles.LOG,'- No file loaded'];
            set(handles.tag_Log,'string',handles.LOG)
            return
        end
        %
        handles.LOG = [handles.LOG,['- Loading ',filename,],[' ... please wait ...']];
        set(handles.tag_Log,'string',handles.LOG)
        %
        [Data,Virtuals_geo]=load_data_from_ensite_mo_v2(filename);
        
        handles.LOG = [handles.LOG,['- Loaded: ',fnamesig]];
        set(handles.tag_Log,'string',handles.LOG)
        
        
        st = {'ECG','EGMs','Virtual'};
        v = listdlg('ListString',st,'SelectionMode','single','PromptString','Select a geometry');
        handles.geoname = st{v};
        %
        handles.LOG = [handles.LOG,['- selected geometry:',st{v}]];
        set(handles.tag_Log,'string',handles.LOG)
        %
        set(handles.tag_name,'string',['- ',fnamesig]);
        ParamSig.Virtuals_geo = Virtuals_geo;
        
        handles.signals_raw=Data.traces;
        ParamSig.frequency = size(Data.traces,1)/(Data.Times(end,1)-Data.Times(1,1));
        ParamSig.Times = Data.Times;
        ParamSig.Times_name = Data.Times_name;
        ParamSig.Type = Data.Type;
        ParamSig.Label = Data.wavenames;
        ParamSig.name_file = fnamesig;
        ParamSig.name_dir = pname;
        ParamSig.geoname = 'Ensite';
        
    case 'Ensite(*.txt)' % Ensite Classic txt
        [fnamesig,pname] = uigetfile('*.txt','Choose a File');
        filename = [pname fnamesig];
        %
        if sum(filename)==0
            set(handles.tag_import_data,'BackgroundColor',[1 1 1]*.941,'fontweight','normal')
            handles.LOG = [handles.LOG,'- No file loaded'];
            set(handles.tag_Log,'string',handles.LOG)
            return
        end
        %
        handles.LOG = [handles.LOG,['- Loading ',filename,],[' ... please wait ...']];
        set(handles.tag_Log,'string',handles.LOG)
        str = {'Waveforms (*Wav*.txt)','Grid of Virtuals (*Virt*.txt)'};
        list_0 = 1;
        if ~isempty(findstr(fnamesig,'Wav'))
            list_0 = 1;
        elseif ~isempty(findstr(fnamesig,'Virt'))
            list_0 = 2 ;
        end
        [s] =listdlg('PromptString','Select type of file:','SelectionMode','single','ListString',str,'InitialValue', list_0);
        
        if s==1
            [signals,ParamOri]=load_data_from_ensite_Classic_Review_Wav(filename);
            handles.signals_raw= signals;
            ParamSig = ParamOri;
            ParamSig.name_file = fnamesig;
            ParamSig.name_dir = pname;
            ParamSig.geoname = 'Ensite_wav';
        elseif s==2
            [signals,ParamOri]=load_data_from_ensite_Classic_Review_Virt(filename);
            handles.signals_raw= signals;
            ParamSig = ParamOri;
            ParamSig.name_file = fnamesig;
            ParamSig.name_dir = pname;
            ParamSig.geoname = 'Ensite_grid';
            
            [fname_geo,pname_geo] = uigetfile('*loc*.txt','Select a ENGUIDE/GEOMETRY file');
            if ~isequal(fname_geo,0)
                filename_geo = [pname_geo,fname_geo];
                [Virtuals_geo,Geometry,Param] = load_data_from_ensite_Classic_Review_Loc(filename_geo);
                ParamSig.Virtuals_geo = Virtuals_geo;
                ParamSig.Geometry = Geometry;
            end
        end
        
        handles.LOG = [handles.LOG,['- Loaded: ',fnamesig]];
        set(handles.tag_Log,'string',handles.LOG)
        %
        set(handles.tag_name,'string',['- ',fnamesig]);
        
        
    case 'MEA(*.txt)'
        [fnamesig,pname] = uigetfile('*.txt', 'Choose a Signal File');
        filename = [pname fnamesig];
        %
        if sum(filename)==0
            set(handles.tag_import_data,'BackgroundColor',[1 1 1]*.941,'fontweight','normal')
            handles.LOG = [handles.LOG,'- No file loaded'];
            set(handles.tag_Log,'string',handles.LOG)
            return
        end
        %
        handles.LOG = [handles.LOG,['- Loading ',filename,],[' ... please wait ...']];
        set(handles.tag_Log,'string',handles.LOG)
        %
        [signals_raw,ParamSig] = Load_MEA_Data(filename);
        
        signals_raw(:,1) = []; % canceling time vector
        ParamSig.Label(1) = [];
        ParamSig.geoname = 'MEA';
        ParamSig.name_file = fnamesig;
        ParamSig.name_dir = pname;
        handles.LOG = [handles.LOG,['- Loaded: ',fnamesig]];
        set(handles.tag_Log,'string',handles.LOG)
        handles.signals_raw=signals_raw;
        
    case 'Biosemi(*.bdf)' % Biosemi data
        [fnamesig,pname] = uigetfile('*.bdf', 'Choose a Signal File');
        filename = [pname fnamesig];
        %
        [signals_raw,~,Label,txt,fs,gain,prefiltering,ChanDim] = eeg_read_bdf(filename,'all','n');
        ParamSig.Label = Label;
        ParamSig.geoname = 'Biosemi';
        ParamSig.name_file = fnamesig;
        ParamSig.name_dir = pname;
        handles.LOG = [handles.LOG,['- Loaded: ',fnamesig]];
        set(handles.tag_Log,'string',handles.LOG)
        handles.signals_raw = detrend(signals_raw','constant');
        ParamSig.frequency = fs;
        ParamSig.gain = gain;
        ParamSig.ChanDim = ChanDim;
        ParamSig.prefiltering = prefiltering;
        
        %%
    case 'GE(*.xml)'
        [fnamesig,dir_load] = uigetfile('*.xml', 'Choose a Signal File');
        fid=fopen([dir_load,fnamesig]);
        fseek(fid,0,'bof');
        %% Full disclosure
        seekstr_strp = '<FullDisclosure>';
        vcontrol = 1;
        while vcontrol
            a=fgetl(fid);
            if ~ischar(a), vcontrol =0; end
            if (strfind(a, seekstr_strp)), break, end
        end
        b = fgetl(fid);
        if  isequal(b(1:5),'Error')
            disp(['=> ',fnamesig, 'is corrupted'])
            figure('visible','off')
            set(gcf,'units','centimeters','position',[2 2 15 12],'paperposition',[2 2 15 12],'papersize',[15 12])
            ant = fnamesig(1:end-4);ant(ant=='_')='-';
            annotation('textbox',[0 .97 1 .03],'string',ant,'horizontalalignment','center','linestyle','none','fontsize',14)
            
            print(gcf,[dir_save,fnamesig(1:end-4),'_Full_Disclosure'],'-djpeg')
            close
            
            Strips =[];
            FullDiscosure = [];
            save([dir_save,fnamesig(1:end-4),'_ECG'],'Strips','FullDiscosure')
            clearvars -except dir_load  isbj dir_save SbjNames do_plot_full
            
            return
        end
        if contains(b,'</FullDisclosure>')
            helpdlg(['Sorry No Full Disclosure for ',fnamesig])
            set(handles.tag_import_data,'BackgroundColor',[1 1 1]*.941,'fontweight','normal')
            error(['Sorry No Full Disclosure for ',fnamesig])
            
        end
        amin=fgetl(fid);
        asec=fgetl(fid);
        imin1 = find(amin=='>');imin2 = find(amin=='<');
        isec1 = find(asec=='>');isec2 = find(asec=='<');
        if ~isempty(imin1)
            FullDiscosure.StartTime(1) = str2double(amin(imin1(1)+1:imin2(2)-1));
            FullDiscosure.StartTime(2) = str2double(asec(isec1(1)+1:isec2(2)-1));
        else
            FullDiscosure.StartTime(1:2)=nan;
        end
        FullDiscosure.TimeCell = {[num2str(FullDiscosure.StartTime(1)),'min ',num2str(FullDiscosure.StartTime(2)),'sec']};
        % N of channels
        fgetl(fid);
        aNchan=fgetl(fid);
        i1 = find(aNchan=='>');i2 = find(aNchan=='<');
        FullDiscosure.Number_of_Channels = str2double(aNchan(i1(1)+1:i2(2)-1));
        % Fs
        aFs=fgetl(fid);
        i1 = find(aFs=='>');i2 = find(aFs=='<');
        FullDiscosure.Frequency_sample = str2double(aFs(i1(1)+1:i2(2)-1));
        i1 = find(aFs=='"');
        FullDiscosure.Frequency_sample_Units = aFs(i1(1)+1:i1(2)-1);
        % Resolution
        aRes=fgetl(fid);
        i1 = find(aRes=='>');i2 = find(aRes=='<');
        FullDiscosure.Resolution = str2double(aRes(i1(1)+1:i2(2)-1));
        i1 = find(aRes=='"');
        FullDiscosure.Resolution_Units = aRes(i1(1)+1:i1(2)-1);
        % Leads
        aLeads=fgetl(fid);
        i1 = find(aLeads=='>');i2 = find(aLeads=='<');
        xx = textscan(aLeads(i1(1)+1:i2(2)-1),'%s ','delimiter',',');
        FullDiscosure.Leads_name = xx{1};
        
        seekstr_strp = '<FullDisclosureData>';
        vcontrol = 1;
        while vcontrol
            a=fgetl(fid);
            if ~ischar(a), vcontrol =0; end
            if (strfind(a, seekstr_strp)), break, end
        end
        data_all = textscan(fid,'%f ','delimiter',',');
        L = length(data_all{1});
        if rem(L,FullDiscosure.Number_of_Channels)~=0
            error('Check full disclosure data in the file')
        end
        y = reshape(data_all{1},FullDiscosure.Frequency_sample,FullDiscosure.Number_of_Channels,length(data_all{1})/FullDiscosure.Frequency_sample/FullDiscosure.Number_of_Channels);
        ecg = nan(size(y,1)*size(y,3),FullDiscosure.Number_of_Channels);
        for i = 1:size(y,3)
            hh = [(FullDiscosure.Frequency_sample)*(i-1)+1 : FullDiscosure.Frequency_sample*i];
            ecg(hh,:) = y(:,:,i);
        end
        clear i hh y
        
        FullDiscosure.ecg = ecg;
        clear a* ecg
        
        %%
        ParamSig.Label = FullDiscosure.Leads_name;
        ParamSig.geoname = 'GE';
        ParamSig.name_file = fnamesig;
        ParamSig.name_dir = dir_load;
        handles.LOG = [handles.LOG,['- Loaded: ',fnamesig]];
        set(handles.tag_Log,'string',handles.LOG)
        handles.signals_raw = detrend(FullDiscosure.ecg,'constant');
        ParamSig.frequency = FullDiscosure.Frequency_sample;
        ParamSig.gain = FullDiscosure.Resolution;
        ParamSig.gain_unit = FullDiscosure.Resolution_Units;
        ParamSig.Times = FullDiscosure.StartTime;
        ParamSig.Times_name = FullDiscosure.TimeCell;
        
    case '*.edf' % EDF files
        
        
        [fnamesig,dir_load] = uigetfile('*.edf', 'Choose a Signal File');
        
        [hdr] = edfread([dir_load,fnamesig]);
        
        Str = cell(hdr.ns,1);
        for i =1:length(hdr.label)
            Str(i,:) = {[num2str(i),' - ',hdr.label{i}]};
        end
        [SELECTION] = listdlg('ListString',Str,'PromptString','Select channels to export:');
        
        [hdr, Stot] = edfread([dir_load,fnamesig],'targetSignals',SELECTION);
        
        
        signals_raw = Stot.';
        
        if isfield(hdr,'samples')
            ParamSig.frequency = hdr.samples;
        else
            answ = inputdlg('Enter sampling frequency [Hz]','sampling frequency',1,{'1000'});
            ParamSig.frequency = str2double(answ{1});
        end
        ParamSig.hdr = hdr;
        ParamSig.Labels = hdr.label;
        ParamSig.geoname = 'EDF';
        ParamSig.name_file = fnamesig;
        ParamSig.name_dir = dir_load;
        handles.LOG = [handles.LOG,['- Loaded: ',fnamesig]];
        set(handles.tag_Log,'string',handles.LOG)
        handles.signals_raw=signals_raw;
        
    case 'CardioInsight'
        % == Load signals
        [fnamesig,dir_load] = uigetfile('*.txt', 'Choose a Signal File');
        ParamSig.name_file = fnamesig;
        fileID = fopen([dir_load,fnamesig],'r');
        a = fgetl(fileID);
        N = sum(a==',')+1;
        fseek(fileID,0,-1);
        formatSpec = repmat('%f',[1 N]);
        dataArray = textscan(fileID, formatSpec, 'Delimiter', ',', 'EmptyValue' ,NaN, 'ReturnOnError',false,'collectoutput',1);
        fclose(fileID);
        signals = [dataArray{1}];
        ParamSig.filename_sig = [dir_load,fnamesig];
        clear dataArray
        strload = {['- Loaded: ',fnamesig]};
        % == Load geometry [Ventricles]
        clear filename
        ii = find(dir_load=='\');
        dir_load_main = dir_load(1:ii(end)-1);
        [filename,dir_load] = uigetfile([dir_load_main,'\*.csv'], 'Choose a geometry File [VENTRICLES]');
        if filename~=0
            name_main = filename;
            ii = find(name_main=='.');
            name_main=name_main(1:ii(2)-1);
            
            T = csvread([dir_load,name_main,'.Faces.csv']);
            geo.Mesh_ventricles = T+1;
            T = csvread([dir_load,name_main,'.Vertices.csv']);
            geo.Nodes_ventricles = T*10; % cm -> mm (double check)
            ParamSig.dir_geo = dir_load;
            strload = cat(1,strload,{['- Loaded: ',filename]});
            
        end
        
        % == Load geometry [LAD]
        clear filename ii T
        [filename,dir_load] = uigetfile([dir_load,'\*.csv'], 'Choose a geometry File [LAD]');
        if filename~=0
            name_main = filename;
            ii = find(name_main=='.');
            name_main=name_main(1:ii(2)-1);
            
            T = csvread([dir_load,name_main,'.Faces.csv']);
            geo.Mesh_LAD = T+1;
            T = csvread([dir_load,name_main,'.Vertices.csv']);
            geo.Nodes_LAD = T*10; % cm -> mm (double check)
            strload = cat(1,strload,{['- Loaded: ',filename]});
        end
        
        % == Load geometry [LAD]
        clear filename T ii
        [filename,dir_load] = uigetfile([dir_load,'\*.csv'], 'Choose a geometry File [VALVES]');
        if filename~=0
            name_main = filename;
            ii = find(name_main=='.');
            name_main=name_main(1:ii(2)-1);
            
            T = csvread([dir_load,name_main,'.Faces.csv']);
            geo.Mesh_Valves = T+1;
            T = csvread([dir_load,name_main,'.Vertices.csv']);
            geo.Nodes_Valves = T*10; % cm -> mm (double check)
            strload = cat(1,strload,{['- Loaded: ',filename]});
            
        end
        
        
        % == Load geometry [LAD]
        clear filename T ii
        [filename,dir_load] = uigetfile([dir_load,'\*.csv'], 'Choose a geometry File [AORTA]');
        if filename~=0
            name_main = filename;
            ii = find(name_main=='.');
            name_main=name_main(1:ii(2)-1);
            
            T = csvread([dir_load,name_main,'.Faces.csv']);
            geo.Mesh_Aorta = T+1;
            T = csvread([dir_load,name_main,'.Vertices.csv']);
            geo.Nodes_Aorta = T*10; % cm -> mm (double check)
            strload = cat(1,strload,{['- Loaded: ',filename]});
            
        end
        
        % == Load geometry [ATRIA]
        clear filename T ii
        [filename,dir_load] = uigetfile([dir_load,'\*.csv'], 'Choose a geometry File [ATRIA]');
        if filename~=0
            name_main = filename;
            ii = find(name_main=='.');
            name_main=name_main(1:ii(2)-1);
            
            T = csvread([dir_load,name_main,'.Faces.csv']);
            geo.Mesh_Atria = T+1;
            T = csvread([dir_load,name_main,'.Vertices.csv']);
            geo.Nodes_Atria = T*10; % cm -> mm (double check)
            strload = cat(1,strload,{['- Loaded: ',filename]});
            
        end
        
        ParamSig.frequency = 1000;
        handles.signals = signals;
        handles.geo = geo;
        handles.LOG = [handles.LOG,strload'];
        set(handles.tag_Log,'string',handles.LOG)
        ParamSig.geoname = 'CardioInsight';
        %         keyboard
        
    case 'CPH (*.dasc)'
        [fnamesig,dir_load] = uigetfile('*.dasc', 'Choose a Signal File');
        [signals,ParamSig] = load_dasc([dir_load,fnamesig]);
        handles.LOG = [handles.LOG,['- Loaded: ',fnamesig]];
        set(handles.tag_Log,'string',handles.LOG)
        ParamSig.name_file = fnamesig;
        
        str = {'new_sock4(Stefan)','old_sock4(Stefan)','none'};
        [s] = listdlg('PromptString','Select a file:',...
            'SelectionMode','single',...
            'ListString',str)
        
        if s==1
            V = load('.\GUI_egm_mFiles\Geo_new_sock4_Stefan');
            signals = signals(:,V.ichan);
            ParamSig.geo = V.geo;
            ParamSig.Label = V.Label;
        elseif s==2
            V = load('.\GUI_egm_mFiles\Geo_old_sock4_Stefan');
            signals = signals(:,V.ichan);
            ParamSig.geo = V.geo;
            ParamSig.Label = V.Label;
        end
        
        handles.signals_raw=signals;
        
    case 'CARTO (*Raw Data)'
        
        addpath .\GUI_Load_Carto
        GUI_Load_Carto;
        close(handles.figure1)
        return
        
    case 'Hlamp (*hdf5)'
        
        [fnamesig,p] = uigetfile('*.hdf5');
        fn = [p,fnamesig];
        signals_raw = h5read(fn,'/RawData/Samples').'; % Read signal / column vector
        signals_raw = double(signals_raw/1000); % mV
        ParamAmp = ReadProperties_h5(fn);
        ParamSig.frequency = ParamAmp.SampleRate;
        ind = 0;
        L = cell(1,256);
        for i = 1:64
            ind = ind+1;
            L{ind} = [num2str(ind),'-A',num2str(i)];
        end
        for i = 1:64
            ind = ind+1;
            L{ind} = [num2str(ind),'-B',num2str(i)];
        end
        for i = 1:64
            ind = ind+1;
            L{ind} = [num2str(ind),'-C',num2str(i)];
        end
        for i = 1:64
            ind = ind+1;
            L{ind} = [num2str(ind),'-D',num2str(i)];
        end
        ParamSig.Label = L;
        ParamSig.geoname = 'HDF5';
        ParamSig.name_file = fnamesig;
        ParamSig.name_dir = p;
        handles.LOG = [handles.LOG,['- Loaded: ',fn]];
        set(handles.tag_Log,'string',handles.LOG)
        
        handles.signals_raw=signals_raw;
        
    case 'Precision(zipped DxL)'
        [fnamesig,p] = uigetfile('*.zip');
        fn = [p,fnamesig];
        OUT = Load_case_Precision_fun_for_GUI(fn);
        
        handles.signals = OUT.S;
        ParamSig.Label = OUT.labels;
        ParamSig.frequency = OUT.fs;
        ParamSig.Map_name = OUT.Map_name;
        geo.Cmesh.triangles = OUT.MESH.MESH.Faces;
        geo.Cmesh.vertices = OUT.MESH.MESH.vertices;
        geo.Cmesh.vertices_norm = OUT.MESH.MESH.Normals_number;
        geo.Study_name = OUT.MESH.Study_name;
        geo.xyz = OUT.xyz.';
        geo.xyz_surfP = OUT.xyz_surfP.';
        handles.geo = geo;
        
        handles.LOG = [handles.LOG,['- Loaded: ',fn]];
        set(handles.tag_Log,'string',handles.LOG)
end
%%
set(handles.tag_name,'string',['- ',fnamesig]);


if ~exist('ParamSig','var')
    ff= inputdlg('Please enter the sampling frequency [Hz]','INPUT',1,{'1000'});
    ParamSig.frequency = str2double(ff);
end
if ~isfield(ParamSig,'Label')
    if exist('signals','var')
        for i=1:size(signals,2)
            ParamSig.Label(i) = {num2str(i)};
        end
    end
end
if ~isfield(ParamSig,'Label')
    if exist('signals_raw','var')
        for i=1:size(signals_raw,2)
            ParamSig.Label(i) = {num2str(i)};
        end
    end
end

handles.ParamSig=ParamSig;
guidata(hObject, handles);
set(handles.tag_import_data,'BackgroundColor',[1 1 1]*.941,'fontweight','normal')



% --- Executes on selection change in tag_Log.
function tag_Log_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_Log contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_Log

% --- Executes during object creation, after setting all properties.
function tag_Log_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tag_import_type.
function tag_import_type_Callback(hObject, eventdata, handles)
% hObject    handle to tag_import_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_import_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_import_type


% --- Executes during object creation, after setting all properties.
function tag_import_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_import_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in tag_save_data.
function tag_save_data_Callback(hObject, eventdata, handles)
% hObject    handle to tag_save_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

name_temp = handles.ParamSig.name_file;
name_temp(find(name_temp=='.'):end)=[];
set(handles.tag_save_data,'backgroundcolor',[.6 .6 .6],'foregroundcolor','r','FontWeight','bold')
v = fieldnames(handles);
v(strncmp(v,'tag*',3))=[];v(strncmp(v,'uipan',5))=[];v(strncmp(v,'text',4))=[];
v(strncmp(v,'type',4))=[];
v(strcmp(v,{'figure1'}))=[];v(strcmp(v,{'output'}))=[];v(strncmp(v,'hs',2))=[];
v(strcmp(v,{'Markers'}))=[];v(strcmp(v,{'LOG'}))=[];
v = sort(v);
iv = listdlg('liststring',v,'PromptString','Select variables to save');
[fname,pathn] = uiputfile('*.mat','Save file name',name_temp);
filename = [pathn,fname];
if sum(filename)~=0
    hw = waitbar(0,'Please wait untli saving is finished ...');
    for i = 1:length(iv)
        eval([v{iv(i)} ,'=handles.',(v{iv(i)}),';']);
        if i==1
            save(filename,v{iv(i)});
        else
            waitbar(i/length(iv),hw);
            
            
            if isequal(v{iv(i)},'Resti_FIT_exp')
                if isfield(Resti_FIT_exp,'curve')
                    Resti_FIT_exp = rmfield(Resti_FIT_exp,'curve');
                    Resti_FIT_exp = rmfield(Resti_FIT_exp,'gof');
                end
            end
            if isequal(v{iv(i)},'Resti_FIT_lin') % just to save time
                if isfield(Resti_FIT_lin,'curve')
                    Resti_FIT_lin = rmfield(Resti_FIT_lin,'curve');
                    Resti_FIT_lin = rmfield(Resti_FIT_lin,'gof');
                end
            end
            save(filename,v{iv(i)},'-append');
        end
    end
    close(hw)
    handles.LOG = [handles.LOG,['- Saved: ',filename],'--------------','- Variables:',v','--------------'];
    set(handles.tag_Log,'string',handles.LOG)
else
    handles.LOG = [handles.LOG,'- Saving aborted by the user'];
    set(handles.tag_Log,'string',handles.LOG)
end
set(handles.tag_save_data,'backgroundcolor',[1 1 1]*.941,'foregroundcolor','k','FontWeight','normal')



% --- Executes on selection change in tag_method.
function tag_method_Callback(hObject, eventdata, handles)
% hObject    handle to tag_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_method



function tag_name_Callback(hObject, eventdata, handles)
% hObject    handle to tag_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_name as text
%        str2double(get(hObject,'String')) returns contents of tag_name as a double


% --- Executes during object creation, after setting all properties.
function tag_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in tag_analysis_special.
function tag_analysis_special_Callback(hObject, eventdata, handles)
% hObject    handle to tag_analysis_special (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.tag_analysis_special,'backgroundColor',[.5 .5 .5],'foregroundColor','r')


S = {'ECGI - CardioInsight','Body Surface Potentials','CARTO','Precision','NavX Classic','Sock'};
[is] = listdlg('ListString',S);
SELECTION = S{is};

switch SELECTION
    case 'NavX Classic'
        
        if isfield(handles.ParamSig,'geo'); % data from Stefan Denmark
            handles.geo = handles.ParamSig.geo;
        else
            geo_new = get_Ensite_mesh(handles.ParamSig.Virtuals_geo.xyz);
            handles.geo = geo_new;
        end
        guidata(hObject,handles)
        
        GUI_EGM_map_viewer_CardioInsight_and_Ensite(hObject)
    case 'ECGI - CardioInsight'
        GUI_EGM_map_viewer_CardioInsight_v3(hObject)
    case 'Body Surface Potentials'
        GUI_BSPM_map_viewer(hObject)
    case 'CARTO'
        GUI_CARTO_viewer(hObject)
    case 'Sock'
        GUI_EGM_map_viewer(hObject)
    case 'Precision'
        GUI_Precision_CPH_viewer(hObject)
end
set(handles.tag_analysis_special,'backgroundColor', [1 1 1]*.941,'foregroundColor','k')



% --- Executes on button press in tag_plot.
function tag_plot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hplot = PlotSignals(hObject);
% guidata(hObject);
if isfield(handles,'signals')&isfield(handles,'signals_raw')
    handles = rmfield(handles,'signals_raw');
    guidata(hObject,handles);
end
%
handles.LOG = [handles.LOG,'- Signals_raw removed to save space'];
set(handles.tag_Log,'string',handles.LOG)
%

% --- Executes on button press in tag_spikes_gui.
function tag_spikes_gui_Callback(hObject, eventdata, handles)
% hObject    handle to tag_spikes_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.tag_spikes_gui,'BackgroundColor',[.6 .6 .6],'fontweight','bold','foregroundColor','r')
handles.LOG = [handles.LOG,['- Semi-automatic localization of spikes ...']];
set(handles.tag_Log,'string',handles.LOG)
%
hspikes = Spikes_GUI(hObject);
set(handles.tag_spikes_gui,'BackgroundColor',[1 1 1]*.941,'fontweight','normal','foregroundColor','k')


% --- Executes on button press in tag_markerrs.
function tag_markerrs_Callback(hObject,eventdata,handles)
% hObject    handle to tag_markerrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% vover = {'signals_raw','signals_proc','S','Markers','MarkersC','SNR'};
% for i=1:length(vover)
%     if isfield(handles,vover{i})
%         handles = rmfield(handles,vover{i});
%     end
% end
% guidata(hObject,handles)

if ~isfield(handles,'geoname')
    handles.geoname = 'cath-lab';
end

% hs_Mcorrect = Markers_GUI(hObject);

if ~isfield(handles,'signals')
    handles.signals = handles.signals_raw;
end

data.S = handles.signals;

data.ParamSig = handles.ParamSig;
data.spikes = handles.spikes(:);
if isfield(handles,'signals_proc')
    data.signals_proc = handles.signals_proc;
end
if isfield(handles,'SNR')
    data.SNR = handles.SNR;
end
if isfield(handles,'Markers')
    data.Markers = handles.Markers;
end
if isfield(handles,'MarkersC')
    data.MarkersC = handles.MarkersC;
end
if isfield(handles,'Markers')&~isfield(data,'MarkersC')
    data.MarkersC = handles.Markers;
end
if isfield(handles,'sig_corr')
    data.sig_corr = handles.sig_corr;
end
if isfield(handles,'sig_corr_beats')
    data.sig_corr_beats = handles.sig_corr_beats;
end
if isfield(handles,'signals_proc_AT')
    data.signals_proc_AT = handles.signals_proc_AT;
end
data.hObject_main = hObject;

hs_Mcorrect = Markers_GUI_v3(data);


% --- Executes on button press in tag_relabel.
function tag_relabel_Callback(hObject, eventdata, handles)
% hObject    handle to tag_relabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hReLab = ReLable_GUI(hObject);



% --- Executes on button press in tag_copy_to_workspace.
function tag_copy_to_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to tag_copy_to_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = handles;
v = fieldnames(data);
iv2 = find([strncmp(v,'tag_',4) | strncmp(v,'text',4) | strncmp(v,'uip',3) | strncmp(v,'figu',4)]);
data = rmfield(data,v(iv2));
v = fieldnames(data);
ii = listdlg('ListString',v,'PromptString','Select variables to pass to the workspace');
data = rmfield(data,v(setdiff([1:length(v)],ii)));
assignin('base','DataFromGui',data);
v = fieldnames(data);
handles.LOG = [handles.LOG,['- DataFromGui copied in base workspace:']];
handles.LOG = [handles.LOG,v'];
set(handles.tag_Log,'string',handles.LOG)
