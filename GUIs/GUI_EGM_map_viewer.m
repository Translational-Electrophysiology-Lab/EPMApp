function varargout = GUI_EGM_map_viewer(varargin)
% GUI_EGM_MAP_VIEWER MATLAB code for GUI_EGM_map_viewer.fig
%      GUI_EGM_MAP_VIEWER, by itself, creates a new GUI_EGM_MAP_VIEWER or raises the existing
%      singleton*.
%
%      H = GUI_EGM_MAP_VIEWER returns the handle to a new GUI_EGM_MAP_VIEWER or the handle to
%      the existing singleton*.
%
%      GUI_EGM_MAP_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_EGM_MAP_VIEWER.M with the given input arguments.
%
%      GUI_EGM_MAP_VIEWER('Property','Value',...) creates a new GUI_EGM_MAP_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_EGM_map_viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_EGM_map_viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_EGM_map_viewer

% Last Modified by GUIDE v2.5 29-Jan-2015 18:45:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_EGM_map_viewer_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_EGM_map_viewer_OutputFcn, ...
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


% --- Executes just before GUI_EGM_map_viewer is made visible.
function GUI_EGM_map_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_EGM_map_viewer (see VARARGIN)

% Choose default command line output for GUI_EGM_map_viewer
handles.output = hObject;


axis(handles.tag_mesh,'off')
axis(handles.tag_egm,'off')
handles.cclines.handles = [];
handles.cclines.color = [];
%
handles.hmain = guidata(varargin{1});
% if ~isfield(handles,'Markers')
%     handles.hmain.MarkersC = handles.hmain.Markers;
% end
% if ~isfield(handles,'MarkersC')
%     handles.hmain.MarkersC = handles.hmain.MarkersC;
% end


set(handles.tag_beat_beging,'string','1');
set(handles.tag_beat_end,'string',num2str(size(handles.hmain.MarkersC.dt,1)));
set(handles.tag_slider_beat_b,'value',1)
set(handles.tag_slider_beat_e,'value',size(handles.hmain.MarkersC.dt,1))
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_EGM_map_viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_EGM_map_viewer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in tag_menu.
function tag_menu_Callback(hObject, eventdata, handles)
% hObject    handle to tag_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_menu


% --- Executes during object creation, after setting all properties.
function tag_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_plot_egm.
function tag_plot_egm_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_egm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode(handles.figure1);
info_struct = getCursorInfo(dcm_obj);
pos = info_struct.Position;
[~,im] = min(sum(abs(handles.MESH.Vertices_Original  - repmat(pos,[size(handles.MESH.Vertices_Original,1),1])),2));
ichan = handles.MESH.channel_num(im);
ichan_title = [num2str(handles.MESH.channel_num(im)),'/',handles.MESH.elect_name{im}];

if isnan(ichan)
    h = msgbox('Channel not found');
    pause(1)
    close(h)
    return 
end

if get(handles.tag_do_filter,'value')
    S = handles.hmain.signals_proc(:,ichan);
else
    S = handles.hmain.signals(:,ichan);
end
t = [1 : size(S,1)]/handles.hmain.ParamSig.frequency*1000;

if get(handles.tag_hold_on,'value')==1
    hold(handles.tag_egm,'on'),
else
    hold(handles.tag_egm,'off'),
    handles.cclines.handles = [];
    handles.cclines.color = [];
    % Update handles structure
    guidata(hObject, handles);
end
plot(handles.tag_egm,t,S)
hold(handles.tag_egm,'on'),
set(handles.tag_egm,'box','off')
xlabel(handles.tag_egm,'Time (ms)')

cc = get(handles.tag_egm,'children');
cclines = cc(strcmp(get(cc,'marker'),'none'));
if ~isempty(findobj(handles.tag_egm,'type','line','tag','timeline'))
%     cclines(cclines ==findobj(handles.tag_egm,'type','line','tag','timeline'))=[];
    cclines = setdiff(cclines,findobj(handles.tag_egm,'type','line','tag','timeline'));
end
cmap = lines;

if isempty(handles.cclines.handles)
    handles.cclines.handles = cclines;
    handles.cclines.color = cmap(1,:);
else
    handles.cclines.handles = [handles.cclines.handles setdiff(cclines(:)',handles.cclines.handles)];
    handles.cclines.color = cmap(1:length(handles.cclines.handles),:);
end
guidata(hObject, handles);

for i = 1:length(cclines)
    set(handles.cclines.handles(i),'color',handles.cclines.color(i,:))
end

cc0 = get(handles.figure1,'children');
ccleg = cc0(strcmp(get(cc0,'tag'),'legend'));
if isempty(ccleg)
    strleg = {ichan_title};
else
    strleg = get(ccleg,'string');
    strleg = cat(1,strleg(:),[num2str(handles.MESH.channel_num(im)),'/',handles.MESH.elect_name{im}]);
end
legend(handles.cclines.handles(:),strleg,'box','off')

dtp = handles.hmain.MarkersC.dt(:,ichan);
dtp(isnan(dtp)) = [];
rtp = handles.hmain.MarkersC.rt_Wyatt(:,ichan);
rtp(isnan(rtp)) = [];
plot(handles.tag_egm,dtp,S(round(dtp/1000*handles.hmain.ParamSig.frequency),:),'or','linewidth',2,'markersize',8)
plot(handles.tag_egm,rtp,S(round(rtp/1000*handles.hmain.ParamSig.frequency),:),'xr','linewidth',2,'markersize',9)


ibeat = str2double(get(handles.tag_beat_beging,'string')):str2double(get(handles.tag_beat_end,'string'));
ibeat = ibeat(1);
if ibeat < length(handles.hmain.spikes)
    t2 = [handles.hmain.spikes(ibeat) : 1/handles.hmain.ParamSig.frequency*1000 : handles.hmain.spikes(ibeat + 1)];
else
    t2 = [handles.hmain.spikes(ibeat) : 1/handles.hmain.ParamSig.frequency*1000 : size(handles.hmain.signals,1)/handles.hmain.ParamSig.frequency*1000];
end
t0 = get(handles.tag_slider,'value')*range(t2) + handles.hmain.spikes(ibeat);
plot(handles.tag_egm,[1 1]*t0,get(handles.tag_egm,'ylim'),'--k','tag','timeline');

axis(handles.tag_egm,'tight')
set(handles.tag_egm,'visible','on')
set(handles.tag_egm,'xlim',t0+[-2 2]*1000)
% legend(cclines,['node=',num2str(im)])


% --- Executes on button press in tag_map.
function tag_map_Callback(hObject, eventdata, handles)
% hObject    handle to tag_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% this comes from Chris function plotScalarSurface

cla(handles.tag_mesh)
% cla(handles.tag_egm)
cb = findobj(handles.figure1,'tag','Colorbar')
delete(cb);
cb = findobj(handles.figure1,'tag','legend')
delete(cb);
handles.cclines.handles = [];
handles.cclines.color = [];
% cc = get(handles.figure1,'children')
% delete(cc)
clear cb 
ibeat = str2double(get(handles.tag_beat_beging,'string')):str2double(get(handles.tag_beat_end,'string'));
str = get(handles.tag_menu,'string');
istr = get(handles.tag_menu,'value');
if isequal(str{istr},'Activ. T')
    values = handles.hmain.MarkersC.dt; 
    s=handles.hmain.spikes(1:size(values,1));
    values = values - repmat(s(:),[1 size(values,2)]);
    clear s
    colorLabel = 'AT (ms)';
    set(handles.tag_slider,'value',0,'SliderStep',[0.01  0.1])
end
if isequal(str{istr},'Repol. T')
    values = handles.hmain.MarkersC.rt_Wyatt; 
    values = values - repmat(handles.hmain.spikes(:),[1 size(values,2)]);
    colorLabel = 'RT (ms)';
    set(handles.tag_slider,'value',0,'SliderStep',[0.01  0.1])
end
if isequal(str{istr},'ARI')
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
     values = handles.hmain.MarkersC.rt_Wyatt - handles.hmain.MarkersC.dt; 
    colorLabel = 'ARI (ms)';
end

%%
if isequal(str{istr},'Potential')
    ibeat = ibeat(1);
    set(handles.tag_beat_end,'string',num2str(ibeat(1)));
    t = [handles.hmain.spikes(ibeat) : 1/handles.hmain.ParamSig.frequency*1000 : handles.hmain.spikes(ibeat + 1)]-handles.hmain.spikes(ibeat);
    set(handles.tag_slider,'SliderStep',[1/range(t) 0.1])
    t0 = get(handles.tag_slider,'value')*t(end);
    [~,n] = min(abs(t-t0));
    
    
    if get(handles.tag_button_video,'value')
        ii = [round(handles.hmain.spikes(ibeat)/1000*handles.hmain.ParamSig.frequency) : round(handles.hmain.spikes(ibeat+1)/1000*handles.hmain.ParamSig.frequency)];
        if get(handles.tag_do_filter,'value')
            s = diff(handles.hmain.signals_proc(ii(1:end),:));
        else
            s = diff(handles.hmain.signals(ii(1:end),:));
        end
        s = s./repmat(abs(min(s(10:end-10,:))),[size(s,1) 1]);
        values = s(n,:);
        colorLabel = ['(dV)/dt Norm (',num2str(t0,3),'ms)'];
    else
        if get(handles.tag_do_filter,'value')
            values =  handles.hmain.signals_proc(n,:)./range(handles.hmain.signals_proc,1);
        else
            values =  handles.hmain.signals(n,:)./range(handles.hmain.signals,1);
        end
        colorLabel = ['Volt Norm (',num2str(t0,3),'ms)'];
    end
    Clim = [-1 1];
end

%% load MAP
if ~isfield(handles.hmain,'geoname')|~isfield(handles,'MESH')
    st = {'cath-lab','mo_sock2','mo_sock1','new_sock4','old_sock4','old_sock6'};
    v = listdlg('ListString',st,'SelectionMode','single','PromptString','Select a geometry');
    handles.hmain.geoname = st{v};
    load(['.\GUI_egm_mFiles\Geo_Chann\ALLgeoDATA_',handles.hmain.geoname,'_mesh.mat'])
    load(['.\GUI_egm_mFiles\Geo_Chann\ALLgeoDATA_',handles.hmain.geoname,'.mat'])
    handles.MESH = MESH;
    handles.MESH.channel_num = channel_num;
    handles.MESH.elect_name = elect_name;
    guidata(hObject, handles);
end
%%
plotFigure = handles.figure1;
plotAxes = handles.tag_mesh;
figure(plotFigure);
dcm_obj = datacursormode(plotFigure);
cla(plotAxes)
axes(plotAxes);
MESH = handles.MESH;
%% Avarage multiple beats or just take one into account
if ~isequal(str{istr},'Potential')
values = nanmean(values(ibeat,:),1);
end
iiko = handles.hmain.SNR<str2double(get(handles.tag_snr_min,'string')) | handles.hmain.SNR>str2double(get(handles.tag_snr_max,'string'));
values(iiko) = nan;
%%
if ~isequal(str{istr},'Potential')
    Clim = [min(values) prctile(values,98)];
end

if get(handles.tag_do_interp,'value')
    do_interp = 1;
else
    do_interp = 0;
end

% [H,index_val_interp] = patch_index_mo([1:length(values)],values(:),'mo_sock2',do_interp,0);
[H,index_val_interp] = patch_index_mo([1:length(values)],values(:),handles.hmain.geoname,do_interp,0);


handles.values_interp = index_val_interp;
handles.values = values;
guidata(hObject, handles);
% Define some default axes properties so that the user doesn't need to specify them every time
x = MESH.Vertices_Original(:,1);
y = MESH.Vertices_Original(:,2);
z = MESH.Vertices_Original(:,3);

scaleFactor = 0.15;
dx = max(MESH.Vertices_Original(:,1))-min(MESH.Vertices_Original(:,1));
xlim = [min(x)-dx*scaleFactor max(x)+dx*scaleFactor];
dy = max(y)-min(y);
ylim = [min(y)-dy*scaleFactor max(y)+dy*scaleFactor];
dz = max(z)-min(z);
zlim = [min(z)-dz*scaleFactor max(z)+dz*scaleFactor];
set(plotAxes,'xlim',xlim,'ylim',ylim,'zlim',zlim,'alimmode','manual','dataaspectratio',[1 1 1]);


set(H.cbar,'Location','southoutside')
title(H.cbar,colorLabel,'fontsize',12);
set(H.cbar,'position',[.02 .75 .4 .025],'xlim',Clim)
rotate3d(plotAxes);

if exist(['\colormap_carto.mat'],'file')>0
    load('colormap_carto')
    colormap(cmap)
else
    cmap = jet;
    cmap = cmap(end:-1:1,:);
    colormap(cmap)
end

if ~get(handles.tag_show_electrodes,'value')
    set(H.cross,'visible','off')
    set(H.points,'visible','off')
end
if ~get(handles.tag_do_interp,'value')
    cmap(1,:) = [.7 .7 .7];
    colormap(cmap)
end
set(handles.tag_mesh,'clim',Clim)
axis(handles.tag_mesh,'off')

oo = findobj(handles.tag_egm,'tag','timeline');
if ~isempty(oo)
    delete(oo)
    plot(handles.tag_egm,[1 1]*handles.hmain.spikes(ibeat),get(handles.tag_egm,'ylim'),'--k','tag','timeline');
%     set(handles.tag_egm,'xlim',[handles.hmain.spikes(ibeat) handles.hmain.spikes(ibeat)+1000],'box','off')
end

% --- Executes on button press in tag_hold_on.
function tag_hold_on_Callback(hObject, eventdata, handles)
% hObject    handle to tag_hold_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_hold_on


% --- Executes on slider movement.
function tag_slider_Callback(hObject, eventdata, handles)
% hObject    handle to tag_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
str = get(handles.tag_menu,'string');
istr = get(handles.tag_menu,'value');
cbar = findobj(handles.figure1,'tag','Colorbar');

if ~isequal(str{istr},'Potential')
    
    ibeat = str2double(get(handles.tag_beat_beging,'string')):str2double(get(handles.tag_beat_end,'string'));
    ibeat = ibeat(1);
    %     t = [handles.hmain.spikes(ibeat) : 1/handles.hmain.ParamSig.frequency*1000 : handles.hmain.spikes(ibeat + 1)]-handles.hmain.spikes(ibeat);
    
    if get(handles.tag_do_interp,'value')
        values = handles.values_interp;
    else
        values = handles.values;
    end
    Clim = [min(values) prctile(values,98)];

    cla(cbar);
    if size(size(colormap),2)==2
        if exist('colormap_carto.mat','file')
            load('colormap_carto.mat')
        else
            cmap = jet;
            cmap(end:-1:1,:) = cmap;
        end
        colormap(cmap)
    end
    v =get(handles.tag_slider,'value');
    t0 = min(values)+range(values)*v + handles.hmain.spikes(ibeat);
  
    if get(handles.tag_button_video,'value')==0
        if v>0
            set(handles.tag_mesh,'clim',[Clim(1) Clim(1)+(Clim(2)-Clim(1))*v])
        end
    else
        po = findobj(handles.figure1,'type','patch');
        po = findobj(po,'facecolor','interp');
        DeltaT = 20;

        v2 = values;
        v2([ values <= ((min(values)+range(values)*v)-DeltaT/2) | values > ((min(values)+range(values)*v)+DeltaT/2)]) = 0;
        set(po,'FaceVertexCData',v2)
        
        if exist('colormap_carto.mat','file')>1
            load('colormap_carto.mat')
            cmap(1,:) = [.6 .6 .6];
        else
            cmap = jet;
            cmap(end:-1:1,:) = cmap;
            cmap(1,:) = [.6 .6 .6]
        end
        colormap(cmap)
        hold(cbar,'on')
        plot(cbar,[1 1]*(min(values)+range(values)*v),get(cbar,'ylim'),'-k','linewidth',2)
        title(cbar,[num2str(min(values)+range(values)*v,3),' +/- ',num2str(DeltaT/2),' ms'])
    end
    set(cbar,'xlim',Clim)
    
    if ~get(handles.tag_show_electrodes,'value')
        po = findobj(handles.tag_mesh,'marker','.');
        set(po,'visible','off')
        po = findobj(handles.tag_mesh,'marker','x');
        set(po,'visible','off')
    else
        po = findobj(handles.tag_mesh,'marker','.');
        set(po,'visible','on')
        po = findobj(handles.tag_mesh,'marker','x');
        set(po,'visible','on')
    end
    
    oo = findobj(handles.tag_egm,'tag','timeline');
    if ~isempty(oo)
        delete(oo)
    end
    plot(handles.tag_egm,[1 1]*t0,get(handles.tag_egm,'ylim'),'--k','tag','timeline');
    set(handles.tag_egm,'xlim',[handles.hmain.spikes(ibeat) handles.hmain.spikes(ibeat)+1000],'box','off')
    
else
    %% POTENTIAL
    ibeat = str2double(get(handles.tag_beat_beging,'string')):str2double(get(handles.tag_beat_end,'string'));
    
    ibeat = ibeat(1);
    set(handles.tag_beat_end,'string',num2str(ibeat(1)));
%     t = [handles.hmain.spikes(ibeat) : 1/handles.hmain.ParamSig.frequency*1000 : handles.hmain.spikes(ibeat + 1)]-handles.hmain.spikes(ibeat);
    t = [handles.hmain.spikes(ibeat) : 1/handles.hmain.ParamSig.frequency*1000 : handles.hmain.spikes(ibeat + 1)];
    set(handles.tag_slider,'SliderStep',[1/range(t) 10/range(t)])
%     t0 = get(handles.tag_slider,'value')*t(end);% + handles.hmain.spikes(ibeat);
    t0 = get(handles.tag_slider,'value')*range(t) + handles.hmain.spikes(ibeat);
    [~,n] = min(abs(t-t0));
    
    if get(handles.tag_button_video,'value')
        ii = [round(handles.hmain.spikes(ibeat)/1000*handles.hmain.ParamSig.frequency) : round(handles.hmain.spikes(ibeat+1)/1000*handles.hmain.ParamSig.frequency)];
        if get(handles.tag_do_filter,'value')
            s = diff(handles.hmain.signals_proc(ii(10:end-10),:));
        else
            s = diff(handles.hmain.signals(ii(10:end-10),:));
        end
        s = s./repmat(abs(min(s)),[size(s,1) 1]);
        values = s(n,:);
        title(cbar,['(dV)/dt Norm (',num2str(t0-handles.hmain.spikes(ibeat),3),'ms)'],'fontsize',12)
    else
        if get(handles.tag_do_filter,'value')
            values =  handles.hmain.signals_proc(n,:)./range(handles.hmain.signals_proc,1);
        else
            values =  handles.hmain.signals(n,:)./range(handles.hmain.signals,1);
        end
        title(cbar,['Volt Norm (',num2str(t0-handles.hmain.spikes(ibeat),3),'ms)'],'fontsize',12)
    end
    
    po = findobj(handles.figure1,'type','patch');
    po = findobj(po,'facecolor','interp');
    % Re-interpolation
    v = nan(size(handles.MESH.Vertices_Original,1),1);
    [~,ia,ib] = intersect(handles.MESH.channel_num,[1:240]);
    v(ia) =  values(ib);
    iiok = ~isnan(v);
    F = scatteredInterpolant(handles.MESH.Vertices_Original(iiok,1),handles.MESH.Vertices_Original(iiok,2),handles.MESH.Vertices_Original(iiok,3),v(iiok));
    data_interp = F([handles.MESH.Vertices_Interp(:,1),handles.MESH.Vertices_Interp(:,2),handles.MESH.Vertices_Interp(:,3)]);
    set(po,'FaceVertexCData',data_interp)
    
    oo = findobj(handles.tag_egm,'tag','timeline');
    if ~isempty(oo)
        delete(oo)
    end
    plot(handles.tag_egm,[1 1]*t0,get(handles.tag_egm,'ylim'),'--k','tag','timeline');
%     set(handles.tag_egm,'xlim',[t(1) t(end)],'box','off')
end
% end


    
   
    

    



% --- Executes during object creation, after setting all properties.
function tag_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function tag_edit_time_Callback(hObject, eventdata, handles)
% hObject    handle to tag_edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_edit_time as text
%        str2double(get(hObject,'String')) returns contents of tag_edit_time as a double

% --- Executes during object creation, after setting all properties.
function tag_edit_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_filename_Callback(hObject, eventdata, handles)
% hObject    handle to tag_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_filename as text
%        str2double(get(hObject,'String')) returns contents of tag_filename as a double


% --- Executes during object creation, after setting all properties.
function tag_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_button_video.
function tag_button_video_Callback(hObject, eventdata, handles)
% hObject    handle to tag_button_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_button_video


% --- Executes on button press in tag_do_interp.
function tag_do_interp_Callback(hObject, eventdata, handles)
% hObject    handle to tag_do_interp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_do_interp


% --- Executes on button press in tag_show_electrodes.
function tag_show_electrodes_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_electrodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_show_electrodes


% --- Executes on slider movement.
function tag_slider_beat_b_Callback(hObject, eventdata, handles)
% hObject    handle to tag_slider_beat_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

nv = get(handles.tag_slider_beat_b,'value');
if nv > size(handles.hmain.MarkersC.dt,1) 
    set(handles.tag_slider_beat_b,'value',str2double(get(handles.tag_beat_beging,'string')));
end
if nv < 0 
    set(handles.tag_slider_beat_b,'value',0);
end
nv = get(handles.tag_slider_beat_b,'value');
if nv <= size(handles.hmain.MarkersC.dt,1) & nv >0 
    set(handles.tag_beat_beging,'string',num2str(nv))
end


% --- Executes during object creation, after setting all properties.
function tag_slider_beat_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_slider_beat_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function tag_beat_beging_Callback(hObject, eventdata, handles)
% hObject    handle to tag_beat_beging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_beat_beging as text
%        str2double(get(hObject,'String')) returns contents of tag_beat_beging as a double


% --- Executes during object creation, after setting all properties.
function tag_beat_beging_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_beat_beging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function tag_slider_beat_e_Callback(hObject, eventdata, handles)
% hObject    handle to tag_slider_beat_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% d = get(handles.tag_slider_beat_e,'sliderStep');
% nv = str2double(get(handles.tag_beat_end,'string')) + d(1);
% if nv > 0 
%     set(handles.tag_beat_end,'string',num2str(nv))
% end
% nv = get(handles.tag_slider_beat_e,'value');

nv = get(handles.tag_slider_beat_e,'value');
if nv < 0 
    set(handles.tag_slider_beat_e,'value',0);
end
if nv > size(handles.hmain.MarkersC.dt,1) 
    set(handles.tag_slider_beat_e,'value',str2double(get(handles.tag_beat_end,'string')));
end

if nv <= size(handles.hmain.MarkersC.dt,1) & nv >0  
    set(handles.tag_beat_end,'string',num2str(nv))
end


% --- Executes during object creation, after setting all properties.
function tag_slider_beat_e_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_slider_beat_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function tag_beat_end_Callback(hObject, eventdata, handles)
% hObject    handle to tag_beat_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_beat_end as text
%        str2double(get(hObject,'String')) returns contents of tag_beat_end as a double


% --- Executes during object creation, after setting all properties.
function tag_beat_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_beat_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_do_filter.
function tag_do_filter_Callback(hObject, eventdata, handles)
% hObject    handle to tag_do_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_do_filter



function tag_snr_min_Callback(hObject, eventdata, handles)
% hObject    handle to tag_snr_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_snr_min as text
%        str2double(get(hObject,'String')) returns contents of tag_snr_min as a double


% --- Executes during object creation, after setting all properties.
function tag_snr_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_snr_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_snr_max_Callback(hObject, eventdata, handles)
% hObject    handle to tag_snr_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_snr_max as text
%        str2double(get(hObject,'String')) returns contents of tag_snr_max as a double


% --- Executes during object creation, after setting all properties.
function tag_snr_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_snr_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
