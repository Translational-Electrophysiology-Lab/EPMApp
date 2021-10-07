function varargout = GUI_EGM_map_viewer_CardioInsight_and_Ensite(varargin)
% GUI_EGM_MAP_VIEWER_CARDIOINSIGHT_AND_ENSITE MATLAB code for GUI_EGM_map_viewer_CardioInsight_and_Ensite.fig
%      GUI_EGM_MAP_VIEWER_CARDIOINSIGHT_AND_ENSITE, by itself, creates a new GUI_EGM_MAP_VIEWER_CARDIOINSIGHT_AND_ENSITE or raises the existing
%      singleton*.
%
%      H = GUI_EGM_MAP_VIEWER_CARDIOINSIGHT_AND_ENSITE returns the handle to a new GUI_EGM_MAP_VIEWER_CARDIOINSIGHT_AND_ENSITE or the handle to
%      the existing singleton*.
%
%      GUI_EGM_MAP_VIEWER_CARDIOINSIGHT_AND_ENSITE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_EGM_MAP_VIEWER_CARDIOINSIGHT_AND_ENSITE.M with the given input arguments.
%
%      GUI_EGM_MAP_VIEWER_CARDIOINSIGHT_AND_ENSITE('Property','Value',...) creates a new GUI_EGM_MAP_VIEWER_CARDIOINSIGHT_AND_ENSITE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_EGM_map_viewer_CardioInsight_and_Ensite_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_EGM_map_viewer_CardioInsight_and_Ensite_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_EGM_map_viewer_CardioInsight_and_Ensite

% Last Modified by GUIDE v2.5 13-Nov-2018 13:20:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_EGM_map_viewer_CardioInsight_and_Ensite_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_EGM_map_viewer_CardioInsight_and_Ensite_OutputFcn, ...
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


% --- Executes just before GUI_EGM_map_viewer_CardioInsight_and_Ensite is made visible.
function GUI_EGM_map_viewer_CardioInsight_and_Ensite_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_EGM_map_viewer_CardioInsight_and_Ensite (see VARARGIN)

% Choose default command line output for GUI_EGM_map_viewer_CardioInsight_and_Ensite
handles.output = hObject;

%%
h = guidata(varargin{1});
handles.matfile.signals_proc = h.signals_proc;
handles.matfile.signals = h.signals;
handles.matfile.spikes = h.spikes;
if isfield(h,'signals_proc_AT')
    handles.matfile.signals_proc_AT = h.signals_proc_AT;
    handles.tag_plot_filt_AT.Value = 1;
else
    handles.tag_plot_filt_AT.Visible = 'off';
    handles.tag_plot_filt_RT.Value = 1;
end
% =
if ~isfield(h,'SNR');
    P = abs(fft(detrend(handles.matfile.signals,'constant'))).^2;
    ff= [1:size(P,1)]/size(P,1)*h.ParamSig.frequency;
    handles.matfile.SNR = 10*log10( sum(P(ff>2&ff<=40,:))./sum(P(ff>40&ff<=100,:)));
    clear P ff
else
    handles.matfile.SNR = h.SNR;
end
% =
if isfield(h,'MarkersC')
    handles.matfile.Markers = h.MarkersC;
else
    handles.matfile.Markers = h.Markers;
end
handles.matfile.ParamSig = h.ParamSig;
% =
if isfield(h,'RVI')
    handles.matfile.RVI = h.RVI;
end
% =
handles.matfile.geo = h.geo;
handles.hObject_main = varargin{1};



axis(handles.tag_mesh,'off')
axis(handles.tag_egm,'off')
handles.cclines.handles = [];
handles.cclines.color = [];
set(handles.tag_find_chan,'string',handles.matfile.ParamSig.Label,'max',length(handles.matfile.ParamSig.Label));
handles.tag_slider_beat.Max = length(handles.matfile.spikes);
handles.tag_slider_beat.Min = 1;
handles.tag_slider_beat.Value = 1;
if length(handles.matfile.spikes)>1
    handles.tag_slider_beat.SliderStep = [1/(length(handles.matfile.spikes)-1) 4/(length(handles.matfile.spikes)-1)];
else
    handles.tag_slider_beat.SliderStep = [0 0];
end

if isfield(h.ParamSig,'name_file')
    handles.tag_filename.String = h.ParamSig.name_file;
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_EGM_map_viewer_CardioInsight_and_Ensite wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_EGM_map_viewer_CardioInsight_and_Ensite_OutputFcn(hObject, eventdata, handles)
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
xyz = handles.matfile.geo.Nodes_ventricles;

cc = findobj(handles.tag_egm,'linestyle','--','color','k');
delete(cc);
dcm_obj = datacursormode(handles.figure1);
info_struct = getCursorInfo(dcm_obj);
if ~isempty(info_struct)
    pos = info_struct.Position;
else
    ic = handles.tag_find_chan.Value;
    pos = xyz(ic,:);
end

d = sqrt(sum( (xyz  - repmat(pos,[size(xyz ,1),1])).^2 ,2));
[~,im] = sort(d,'ascend');

im = im(1:str2double(get(handles.tag_plot_egm_numb_chan,'string')));

% Do not show those that have been eliminated
if handles.tag_include_NaNs_plot.Value==0
    iiko =  mean(isnan(handles.matfile.Markers.rt_Wyatt(:,im))&isnan(handles.matfile.Markers.dt(:,im)))==1; % 13/11/2018
    im(iiko) = [];
end
if isempty(im);
    cla(handles.tag_egm);
    return
end
set(handles.tag_find_chan,'value',im(1));
handles.tag_r_plot.String = num2str(max(d(im)),'%1.1f');



t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);

%
if handles.tag_plot_filt_RT.Value
    S = handles.matfile.signals_proc(:,im);
    %     handles.tag_plot_filt_RT.Value = 0;
elseif handles.tag_plot_filt_AT.Value
    S = handles.matfile.signals_proc_AT(:,im);
else
    S = handles.matfile.signals(:,im);
    %     handles.tag_plot_filt_RT.Value = 1;
end

ib = str2double(get(handles.tag_beat_number,'string'));

if get(handles.tag_hold_on,'value')==1
    hold(handles.tag_egm,'on'),
    bb = findobj(handles.tag_egm,'tag','UEG');
    xx_lim = get(handles.tag_egm,'xlim');
else
    xx_lim = [handles.matfile.spikes(ib)-100 handles.matfile.spikes(ib)+1000]-handles.matfile.spikes(1);
    cla(handles.tag_egm)
    hold(handles.tag_egm,'off'),
    handles.cclines.handles = [];
    handles.cclines.color = [];
    bb = findobj(handles.tag_egm,'tag','UEG');
    if ~isempty(bb)
        delete(bb);
        bb = [];
    end
    % Update handles structure
    guidata(hObject, handles);
end
bbl = {};
if ~isempty(bb)
    for j = 1:length(bb)
        bbl(j) = {bb(j).DisplayName};
    end
end

aueg = plot(handles.tag_egm,t,S);
set(handles.tag_egm,'xlim',[t(1) t(end)])
zoom(handles.tag_egm,'reset')


% legend(aueg,strleg)
hold(handles.tag_egm,'on'),
set(handles.tag_egm,'box','off')
xlabel(handles.tag_egm,'Time (ms)')

dt = handles.matfile.Markers.dt(ib,im) - handles.matfile.spikes(1);
rt = handles.matfile.Markers.rt_Wyatt(ib,im) - handles.matfile.spikes(1);
for j = 1:length(im)
    % AT
    if ~isnan(dt(j))
        [~,dtp] = min(abs(t-dt(j)));
        plot(handles.tag_egm,t(dtp),S(dtp,j),'marker','o','markerfacecolor',[1 1 1]*.6,'markeredgecolor','k')
    end
    % RT
    if ~isnan(rt(j))
        [~,rtp] = min(abs(t-rt(j)));
        plot(handles.tag_egm,t(rtp),S(rtp,j),'marker','x','markerfacecolor',[1 1 1]*.6,'markeredgecolor','k','linewidth',2)
    end
end
axis(handles.tag_egm,'tight')
set(handles.tag_egm,'visible','on')

if ~get(handles.tag_hold_on,'value')
    aa = findobj(handles.tag_mesh,'marker','x','markersize',15);
    delete(aa);clear aa
end
hold(handles.tag_mesh,'on');

hold(handles.tag_mesh,'on')
for j = 1:length(im)
    aa(j) = plot3(handles.tag_mesh,xyz(im(j),1),xyz(im(j),2),xyz(im(j),3),'x');
    set(aa(j),'markerfacecolor',[1 1 1]*.6,'markeredgecolor','w','markersize',15,'linewidth',2)
end

if ~handles.tag_hold_on.Value | isempty(findobj(handles.tag_egm,'type','line','linestyle','--','color','k'))
    for i= 1:length(handles.matfile.spikes)
        plot(handles.tag_egm,[handles.matfile.spikes(i)]*[1 1]-handles.matfile.spikes(1),get(handles.tag_egm,'ylim'),'--k')
    end
end

plot(handles.tag_egm,[handles.matfile.spikes(ib)]*[1 1]-handles.matfile.spikes(1),get(handles.tag_egm,'ylim'),'--r','linewidth',2)

if ~isequal(xx_lim,[0 1])
    set(handles.tag_egm,'xlim',xx_lim);
    zoom off
end
strleg = cell(size(aueg));
for i= 1:length(aueg)
    strleg(i) = {['node=',num2str(im(i))]};
end
legend([aueg;bb],[strleg(:);bbl(:)])
set(aueg,'tag','UEG')

% cc = get(handles.tag_egm,'children');
% cclines = cc(strcmp(get(cc,'marker'),'none'));
% if ~isempty(findobj(handles.tag_egm,'type','line','tag','timeline'))
%     cclines(cclines ==findobj(handles.tag_egm,'type','line','tag','timeline'))=[];
% end
% if ~isempty(findobj(handles.tag_egm,'type','line','tag','example'))
%     cclines(cclines ==findobj(handles.tag_egm,'type','line','tag','example'))=[];
% end
%
% if length(im)<=64
% cmap = lines;
% else
% cmap = lines(length(im));
% end
% if isempty(handles.cclines.handles)
%     handles.cclines.handles = cclines;
%     handles.cclines.color = cmap(1:length(im),:);
% else
%     handles.cclines.handles = [handles.cclines.handles;setdiff(cclines(:)',handles.cclines.handles)'];
%     handles.cclines.color = cmap(1:length(handles.cclines.handles),:);
% end
% guidata(hObject, handles);
%
% for i = 1:length(cclines)
%     set(handles.cclines.handles(i),'color',handles.cclines.color(i,:))
% end
%
% cc0 = get(handles.figure1,'children');
% ccleg = cc0(strcmp(get(cc0,'tag'),'legend'));
% if isempty(ccleg)
%     for i = 1:length(im)
%         strleg(i) = {['node=',num2str(im(i))]};
%     end
% else
%     strleg = get(ccleg,'string');
%     for i = 1:length(im)
%     strleg = cat(1,strleg(:),['node=',num2str(im(i))]);
%     end
% end
% legend(handles.cclines.handles(:),strleg,'box','off')

% dt = handles.matfile.Markers.dt(im) - handles.matfile.spikes(1);
% rt = handles.matfile.Markers.rt_Wyatt(im) - handles.matfile.spikes(1);
% for j = 1:length(im)
% % AT
% [~,dtp] = min(abs(t-dt(j)));
% plot(handles.tag_egm,t(dtp),S(dtp,j),'marker','o','markerfacecolor',handles.cclines.color(end-length(im)+j,:),'markeredgecolor','k')
% % RT
% [~,rtp] = min(abs(t-rt(j)));
% plot(handles.tag_egm,t(rtp),S(rtp,j),'marker','o','markerfacecolor',handles.cclines.color(end-length(im)+j,:),'markeredgecolor','k')
% end
% axis(handles.tag_egm,'tight')
% set(handles.tag_egm,'visible','on')

% --
% if ~get(handles.tag_hold_on,'value')
%     aa = findobj(handles.tag_mesh,'marker','o');
%     delete(aa);clear aa
% end
%
% for j = 1:length(im)
% aa = plot3(handles.matfile.geo.Nodes_ventricles(im(j),1),handles.matfile.geo.Nodes_ventricles(im(j),2),handles.matfile.geo.Nodes_ventricles(im(j),3),'o');
% set(aa,'markerfacecolor',handles.cclines.color(end-length(im)+j,:),'markeredgecolor','w','markersize',13,'linewidth',2)
% end

% legend(cclines,['node=',num2str(im)])

% --- Executes on button press in tag_load.
function tag_load_Callback(hObject, eventdata, handles)
% hObject    handle to tag_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fnamesig,pname] = uigetfile('*.mat', 'Choose a Signal File');
filename = [pname fnamesig];
V = load(filename);
name = fieldnames(V);
handles.matfile = V;
set(handles.tag_filename,'string',fnamesig,'backgroundcolor',[1 1 1])
%
guidata(hObject,handles);

% --- Executes on button press in tag_map.
function tag_map_Callback(hObject, eventdata, handles)
% hObject    handle to tag_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.tag_interpol.Value=0;
delete(findobj(handles.figure1,'tag','legend'))
delete(findobj(handles.tag_egm,'type','line'))
handles.cclines.handles = [];
handles.cclines.color = [];
cla(handles.tag_mesh)
cla(handles.tag_egm)
guidata(hObject,handles)
set(handles.tag_alpha,'value',0)

ib = str2double(get(handles.tag_beat_number,'string'));
str = get(handles.tag_menu,'string');
istr = get(handles.tag_menu,'value');
if isequal(str{istr},'Activ. T')
    values = handles.matfile.Markers.dt(ib,:);
    values = values - handles.matfile.spikes(ib);
    colorLabel = 'AT (ms)';
    set(handles.tag_slider,'value',0,'SliderStep',[0.01  0.1])
    handles.tag_modify_type.Value = 1;
end
if isequal(str{istr},'Repol. T')
    values = handles.matfile.Markers.rt_Wyatt(ib,:);
    values = values -  handles.matfile.spikes(ib);
    colorLabel = 'RT (ms)';
    set(handles.tag_slider,'value',0,'SliderStep',[0.01  0.1])
    handles.tag_modify_type.Value = 2;
end
if isequal(str{istr},'ARI')
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Markers.rt_Wyatt(ib,:) -handles.matfile.Markers.dt(ib,:);
    colorLabel = 'ARI (ms)';
end



if isequal(str{istr},'Potential')
    t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);
    set(handles.tag_slider,'SliderStep',[1/range(t) 0.1])
    t0 = get(handles.tag_slider,'value')*t(end);
    [~,n] = min(abs(t-t0));
    values =  handles.matfile.signals(n,:)/max(abs(handles.matfile.signals(n,:)));
    colorLabel = ['Volt Norm (',num2str(t0,'%1.0f'),'ms)'];
end
if isequal(str{istr},'dV/dt')
    t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);
    set(handles.tag_slider,'SliderStep',[1/range(t) 0.1])
    t0 = get(handles.tag_slider,'value')*t(end);
    [~,n] = min(abs(t-t0));
    
    s = diff(handles.matfile.signals);
    s = s./repmat(abs(min(s)),[size(s,1) 1]);
    values = s(n,:);
    colorLabel = ['dV/dt Norm (',num2str(t0,'%1.0f'),'ms)'];
    
end

if isequal(str{istr},'RVI')
    D = handles.tag_RVI_d.String;
    rvi_type = handles.tag_RVI_type.Value;
    if rvi_type==1
        name = ['min_D_',D];
    elseif rvi_type==2
        name = ['avg_D_',D];
    end
    
    if isfield( handles.matfile.RVI,name)
        values = handles.matfile.RVI.(name)(ib,:);
    else
        hms = msgbox(['RVI ',name,' has not been calculated']);
        return
    end
    
    if mean(isnan(values))>0.99
        hms = msgbox(['RVI ',name,' has not been calculated for beat = ',num2str(ib)]);
        return
    end
    
    if rvi_type==1
        colorLabel = ['RVI_{min} [D=',num2str(D),'mm]'];
    elseif rvi_type==2
        colorLabel = ['RVI_{avg} [D=',num2str(D),'mm]'];
    end
    
end

if isequal(str{istr},'DispAT')
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_AT(ib,:);
    colorLabel = 'Disp(AT) (ms)';
end
if isequal(str{istr},'DispRT')
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_RT(ib,:);
    colorLabel = 'Disp(RT) (ms)';
end
if isequal(str{istr},'DispARI')
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_ARI(ib,:);
    colorLabel = 'Disp(ARI) (ms)';
end
if isequal(str{istr},'DispAT(N)')
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_AT_norm(ib,:);
    colorLabel = 'Disp(AT) (n.u.)';
end
if isequal(str{istr},'DispRT(N)')
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_RT_norm(ib,:);
    colorLabel = 'Disp(RT) (n.u.)';
end
if isequal(str{istr},'DispARI(N)')
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_ARI_norm(ib,:);
    colorLabel = 'Disp(ARI) (n.u.)';
end


Clim = [prctile(values,5) prctile(values,95)];
if sum(isnan(Clim))>0
    warning('All nan')
    return
end

dcm_obj = datacursormode(handles.figure1);
set(dcm_obj,'UpdateFcn',@myfunctioncursor_markers_EnsiteClassic,'SnapToDataVertex','on')
cla(handles.figure1)
axes(handles.tag_mesh);

%% Meshes 2
MV.faces = handles.matfile.geo.Mesh_ventricles;
MV.vertices =  handles.matfile.geo.Nodes_ventricles;
MV.faceVertexCData = values(:);
%
%% Modify the beginning of this function to just update the colors
hp = findobj(handles.tag_mesh,'type','patch','tag','Mesh_Ventricles');
if ~isempty(hp)
    hp.FaceVertexCData = MV.faceVertexCData;
else
    %% Ventricles
    meshObject = patch(MV,'parent',handles.tag_mesh);
    meshObject.Tag = 'Mesh_Ventricles';
    set(meshObject,'edgecolor','none','facecolor','interp','facelighting','gouraud',...
        'ambientstrength',0.9,'specularcolorreflectance',0.4,'specularexponent',3,'specularstrength',0.3);
    light('Position',[0 0 1],'Style','infinite')
    
    % other structures
    mesh_types = {'Aorta','Atria','PT','septum','LAD','Valves'};
    
    meshObject_other = gobjects(1,length(mesh_types));
    for im = 1:length(mesh_types);
        if isfield(handles.matfile.geo,['Mesh_',mesh_types{im}])
            
            eval(['MV_',mesh_types{im},'.faces = handles.matfile.geo.Mesh_',mesh_types{im},';'])
            eval(['MV_',mesh_types{im},'.vertices = handles.matfile.geo.Nodes_',mesh_types{im},';'])
            eval(['MV_',mesh_types{im},'.faceVertexCData = 1;'])
            
            xx = eval(['MV_',mesh_types{im},';']);
            meshObject_other(im) = patch(xx,'parent',handles.tag_mesh);
            meshObject_other(im).Tag = mesh_types{im};
            set(meshObject_other(im),'edgecolor','k','facecolor',[1 1 1]*im/(1.25*length(mesh_types)));
        else
            eval(['MV_',mesh_types{im},' = [];'])
        end
    end
    
    
    colorbar off
    cbar = findobj(handles.figure1,'tag','Colorbar');
    if isempty(cbar)
        cbar = colorbar;
    end
    
    set(cbar,'Location','southoutside')
    title(cbar,colorLabel,'fontsize',10);
    set(cbar,'position',[.02 .05 .36 .025],'xlim',Clim)
    set(handles.tag_edit_CBmin,'string',num2str(Clim(1),'%1.0f'));
    set(handles.tag_edit_CBmax,'string',num2str(Clim(2),'%1.0f'));
    % camorbit(50,-10);
    rotate3d(handles.tag_mesh);
    
    if exist([pwd,'\GUI_egm_mFiles\Geo_Chann\colormap_carto.mat'],'file')>0
        load([pwd,'\GUI_egm_mFiles\Geo_Chann\colormap_carto.mat'])
        colormap(cmap)
    else
        cmap = jet;
        cmap = cmap(end:-1:1,:);
        colormap(cmap)
    end
    
    set(handles.tag_mesh,'clim',Clim)
    axis(handles.tag_mesh,'off','image')
    
    if get(handles.tag_nodes,'value')
        set(meshObject,'marker','o','MarkerFaceColor',[1 1 1]*.4,'markersize',4)
    else
        set(meshObject,'marker','none','MarkerFaceColor','none')
    end
    
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


if isequal(str{istr},'Activ. T') | isequal(str{istr},'Repol. T') | isequal(str{istr},'ARI')
    %% ACTIVATION TIME
    values = handles.matfile.nodeparam.timestamps.activationTimes;
    values = values - handles.matfile.epipots.t(1);
    Clim = [min(values) prctile(values,98)];
    
elseif isequal(str{istr},'Repol. T')
    %% REPOL TIME
    values = handles.matfile.nodeparam.timestamps.recoveryTimes;
    values = values - handles.matfile.epipots.t(1);
    Clim = [min(values) prctile(values,98)];
elseif isequal(str{istr},'ARI')
    %% ARI
    values = handles.matfile.nodeparam.timestamps.recoveryTimes - handles.matfile.nodeparam.timestamps.activationTimes;
    Clim = [min(values) prctile(values,98)];
end

if isequal(str{istr},'AT')|isequal(str{istr},'ARI')|isequal(str{istr},'RT')
    cla(cbar);
    if size(size(colormap),2)==2
        if exist('colormap_carto.mat','file')
            load('colormap_carto.mat')
        else
            cmap = jet;
            cmap(end:-1:1,:) = cmap;
        end
        po = findobj(handles.figure1,'type','patch');
        po = findobj(po,'facecolor','interp');
        set(po,'FaceVertexCData',values)
        colormap(cmap)
    end
    
    if get(handles.tag_button_video,'value')==0
        
        v =get(handles.tag_slider,'value');
        if v>0
            set(handles.tag_mesh,'clim',[Clim(1) Clim(1)+(Clim(2)-Clim(1))*v])
        end
        
    else
        po = findobj(handles.figure1,'type','patch');
        po = findobj(po,'facecolor','interp');
        DeltaT = 8;
        v =get(handles.tag_slider,'value');
        v2 = values;
        v2([ values <= ((min(values)+range(values)*v)-DeltaT/2) | values > ((min(values)+range(values)*v)+DeltaT/2)]) = 0;
        set(po,'FaceVertexCData',v2)
        
        if exist('colormap_carto.mat','file')
            load('colormap_carto.mat')
            cmap(1,:) = [.6 .6 .6];
        else
            cmap = jet;
            cmap(end:-1:1,:) = cmap;
            cmap(1,:) = [.6 .6 .6]
        end
        colormap(cmap)
        hold(cbar,'on')
        plot(cbar,get(cbar,'xlim'),[1 1]*(min(values)+range(values)*v),'-k','linewidth',2)
        xlabel(cbar,[num2str(min(values)+range(values)*v,3),' +/- ',num2str(DeltaT/2),' ms'])
        
    end
    set(cbar,'xlim',Clim)
    
end

%% POTENTIAL
if isequal(str{istr},'Potential')
    t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);
    set(handles.tag_slider,'SliderStep',[1/range(t) 0.1])
    t0 = get(handles.tag_slider,'value')*t(end);
    [~,n] = min(abs(t-t0));
    values =  handles.matfile.signals(n,:)/prctile(max(abs(handles.matfile.signals)),90);
    colorLabel = ['Volt Norm (',num2str(t0,'%1.0f'),'ms)'];
end

if isequal(str{istr},'dV/dt')
    t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);
    set(handles.tag_slider,'SliderStep',[1/range(t) 0.1])
    t0 = get(handles.tag_slider,'value')*t(end);
    [~,n] = min(abs(t-t0));
    %     s = diff(handles.matfile.signals);
    %     s = s./prctile(abs(min(s)),95);
    %     values = s(n,:);
    s = diff(handles.matfile.signals);
    s = s./repmat(abs(min(s)),[size(s,1) 1]);
    values = s(n,:);
    colorLabel = ['dV/dt Norm (',num2str(t0,'%1.0f'),'ms)'];
    Clim = [prctile(values,5) prctile(values,95)];
end

if isempty(findobj(handles.tag_egm,'type','line','linestyle','-'))
    im = 1;
    S = handles.matfile.signals(:,im);
    plot(handles.tag_egm,t,S,'k');
    aa = plot3(handles.matfile.geo.Nodes_ventricles(im,1),handles.matfile.geo.Nodes_ventricles(im,2),handles.matfile.geo.Nodes_ventricles(im,3),'o');
    set(aa,'markerfacecolor','k','markersize',12)
end

po = findobj(handles.figure1,'type','patch');
po = findobj(po,'facecolor','interp');
set(po,'FaceVertexCData',values(:))

oo = findobj(handles.tag_egm,'tag','timeline');
if ~isempty(oo)
    delete(oo)
end
plot(handles.tag_egm,[1 1]*t0,get(handles.tag_egm,'ylim'),'--k','tag','timeline');
set(handles.tag_egm,'xlim',[t(1) t(end)],'box','off')
title(cbar,colorLabel,'fontsize',10);


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


% --- Executes on button press in tag_modify.
function tag_modify_Callback(hObject, eventdata, handles)
% hObject    handle to tag_modify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cc = findobj(handles.tag_egm,'linestyle','--','color','k');
delete(cc);

dcm_obj = datacursormode(handles.figure1);
info_struct = getCursorInfo(dcm_obj);
pos = info_struct.Position;
samp = round([pos(1)+handles.matfile.spikes(1)]/1000*handles.matfile.ParamSig.frequency);

if get(handles.tag_all_chan,'value')
    aa = findobj(handles.tag_egm,'tag','UEG');
    for i = 1:length(aa)
        ic(i) = str2double(aa(i).DisplayName(6:end));
    end
else
    ic1 = find(handles.matfile.signals(samp,:)==pos(2));
    %     if isfield(handles.matfile,'signals_proc_Ensite_AT')
    if isfield(handles.matfile,'signals_proc_AT')
        ic2 = find(handles.matfile.signals_proc_AT(samp,:)==pos(2));
    else
        ic2 = nan;
    end
    if isfield(handles.matfile,'signals_proc_Ensite_RT')
        ic4 = find(handles.matfile.signals_proc_Ensite_RT(samp,:)==pos(2));
    else
        ic4 = nan;
    end
    ic3 = find(handles.matfile.signals_proc(samp,:)==pos(2));
    
    ic = min([ic1 ic2 ic3 ic4]);
end


if isempty(ic)
    warning('Channel not identified')
    return
end

ib = str2double(handles.tag_beat_number.String);
T = round(str2double(get(handles.tag_window,'string'))/1000*handles.matfile.ParamSig.frequency/2);

% %% modified 29/11/2016
% W = round([pos(1)+handles.matfile.spikes(1)-handles.matfile.spikes(ib)]/1000*handles.matfile.ParamSig.frequency)+[-T:T];
% W(W<1)=1;
% W(W>size(handles.matfile.signals,1))=size(handles.matfile.signals,1);
%
% if ib < length(handles.matfile.spikes);
%     Ws = round(handles.matfile.spikes(ib)/1000*handles.matfile.ParamSig.frequency):round(handles.matfile.spikes(ib+1)/1000*handles.matfile.ParamSig.frequency);
% else
%     Ws = round(handles.matfile.spikes(ib)/1000*handles.matfile.ParamSig.frequency):size(handles.matfile.signals,1);
% end
t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);
Ws = samp + [-T : T]; % modified 29/11/2016
Ws(Ws<1)=[];Ws(Ws>size(handles.matfile.signals,1))=[];
if get(handles.tag_modify_type,'value')==1
    
    if handles.tag_plot_filt_RT.Value
        x = handles.matfile.signals_proc(:,ic);
    elseif handles.tag_plot_filt_AT.Value
        x = handles.matfile.signals_proc_AT(:,ic);
    else
        x = handles.matfile.signals(:,ic);
    end
    X = x(Ws,:);
    
    % Modify activation
    
    Xd = diff(X);Xd2 = diff(Xd);Xd3 = diff(Xd2);
    %     ii = Xd2([W(1):W(end)-1],:).*Xd2([W(1)+1 : W(end)],:)<0 & Xd([W(1): W(end)-1],:)<0 & Xd3([W(1): W(end)-1],:)>0; % zeros of II derivative & I derivative negative & III derivative > 0
    ii = Xd2(1:end-1,:).*Xd2([2 : end],:)<0 & Xd([1: end-2],:)<0 & Xd3([1: end],:)>0; % zeros of II derivative & I derivative negative & III derivative > 0
    
    for j = 1:size(ii,2)
        tt = find(ii(:,j));
        if ~isempty(tt)
            [~,kk]=min(Xd(tt,j));
            %             handles.matfile.Markers.dt(ib,ic(j)) = [tt(kk)+W(1)+Ws(1)-1]/handles.matfile.ParamSig.frequency*1000;
            handles.matfile.Markers.dt(ib,ic(j)) = [tt(kk)+Ws(1)-1]/handles.matfile.ParamSig.frequency*1000;
            % plot new marker
            dt = handles.matfile.Markers.dt(ib,ic(j)) - handles.matfile.spikes(1);
            % AT
            dtp = round(handles.matfile.Markers.dt(ib,ic(j))/1000*handles.matfile.ParamSig.frequency);
            %             plot(handles.tag_egm,t(dtp),handles.matfile.signals(dtp,ic(j)),'xr','markersize',12)
            plot(handles.tag_egm,t(dtp),x(dtp,j),'xr','markersize',12)
        end
    end
    plot(handles.tag_egm,pos(1)-str2double(get(handles.tag_window,'string'))/2*[1 1],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
    plot(handles.tag_egm,pos(1)+str2double(get(handles.tag_window,'string'))/2*[1 1],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
    
    %
    %     po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
    %     cb = get(handles.tag_mesh,'clim');
    %     po.FaceVertexCData = handles.matfile.Markers.dt(ib,:)'-handles.matfile.spikes(ib);
    %     set(handles.tag_mesh,'clim',cb);
    clear x xd*
else
    if handles.tag_plot_filt_RT.Value
        x = handles.matfile.signals_proc(:,ic);
    else
        x = handles.matfile.signals(:,ic);
    end
    
    %     x = handles.matfile.signals_proc(:,ic);
    X = x(Ws,:);
    Xd = diff(X);Xd2 = diff(Xd);Xd3 = diff(Xd2);
    
    %     ii = Xd2([W(1):W(end)-1],:).*Xd2([W(1)+1 : W(end)],:)<0 & Xd([W(1): W(end)-1],:)>0 & Xd3([W(1): W(end)-1],:)<=0;
    ii = Xd2([1:end-1],:).*Xd2([2 : end],:)<0 & Xd([1: end-2],:)>0 & Xd3(1:end,:)<=0;
    
    for j = 1:size(ii,2)
        tt = find(ii(:,j));
        
        if ~isempty(tt)
            [~,kk]=max(Xd(tt,j));
            handles.matfile.Markers.rt_Wyatt(ib,ic(j)) = [tt(kk)+Ws(1)-1]/handles.matfile.ParamSig.frequency*1000;
            
            % plot new marker
            rt = handles.matfile.Markers.rt_Wyatt(ib,ic(j)) - handles.matfile.spikes(1);
            % RT
            %         t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);
            %             [~,dtp] = min(abs(t-rt));
            rtp = round(handles.matfile.Markers.rt_Wyatt(ib,ic(j))/1000*handles.matfile.ParamSig.frequency)+1;
            %             plot(handles.tag_egm,t(dtp),handles.matfile.signals_proc(rtp,ic(j)),'xr','markersize',12)
            plot(handles.tag_egm,t(rtp),x(rtp,j),'xr','markersize',12)
        end
    end
    plot(handles.tag_egm,pos(1)-str2double(get(handles.tag_window,'string'))/2*[1 1],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
    plot(handles.tag_egm,pos(1)+str2double(get(handles.tag_window,'string'))/2*[1 1],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
    
    %
    %     po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
    %     cb = get(handles.tag_mesh,'clim');
    %     po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.spikes(ib);
    %     set(handles.tag_mesh,'clim',cb);
    clear x xd*
end

po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
cb = get(handles.tag_mesh,'clim');
if handles.tag_menu.Value==1
    po.FaceVertexCData = handles.matfile.Markers.dt(ib,:)'-handles.matfile.spikes(ib);
elseif handles.tag_menu.Value==2
    po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.spikes(ib);
elseif handles.tag_menu.Value==3
    po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.Markers.dt(ib,:)';
    
end
set(handles.tag_mesh,'clim',cb);

guidata(hObject,handles)
delete(findall(gcf,'Type','hggroup'));



% --- Executes on selection change in tag_modify_type.
function tag_modify_type_Callback(hObject, eventdata, handles)
% hObject    handle to tag_modify_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_modify_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_modify_type
if handles.tag_modify_type.Value==1
    handles.tag_plot_filt_AT.Value = 1;
    handles.tag_plot_filt_RT.Value = 0;
else
    handles.tag_plot_filt_AT.Value = 0;
    handles.tag_plot_filt_RT.Value = 1;
end

tag_plot_egm_Callback(hObject, eventdata, handles);



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



function tag_window_Callback(hObject, eventdata, handles)
% hObject    handle to tag_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_window as text
%        str2double(get(hObject,'String')) returns contents of tag_window as a double


% --- Executes during object creation, after setting all properties.
function tag_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_nodes.
function tag_nodes_Callback(hObject, eventdata, handles)
% hObject    handle to tag_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_nodes

% cc = findobj(handles.tag_mesh,'type','patch');
% for i = 1:length(cc)
%     L(i) = isequal(cc(i).FaceColor,'interp');
% end
% [~,im] = max(L);
% if get(handles.tag_nodes,'value')
%     set(cc(im),'marker','o','markersize',4,'markerfacecolor','k');
% else
%     set(cc(im),'marker','o','markersize',4,'markerfacecolor','none');
%
% end



% cc = findobj(handles.tag_mesh,'type','line');
%
% if get(handles.tag_nodes,'value')
%     if ~isempty(cc)
%         delete(cc)
%     end
%     hold(handles.tag_mesh,'on');
%     xyz = handles.matfile.geo.Nodes_ventricles;
%     plot3(xyz(:,1),xyz(:,2),xyz(:,3),'ok','parent',handles.tag_mesh,'markerfacecolor',[1 1 1]*.5)
%
% else
%     if ~isempty(cc)
%         delete(cc)
%     end
% end
po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
if get(handles.tag_nodes,'value')
    set(po,'marker','o','MarkerFaceColor',[1 1 1]*.4,'markersize',4)
else
    set(po,'marker','none','MarkerFaceColor','none')
end


% --- Executes on button press in tag_delete_marker.
function tag_delete_marker_Callback(hObject, eventdata, handles)
% hObject    handle to tag_delete_marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode(handles.figure1);

info_struct = getCursorInfo(dcm_obj);
ib = str2double(handles.tag_beat_number.String);

if isequal(get(info_struct.Target.Parent,'tag'),'tag_mesh')
    pos = info_struct.Position;
    xyz = handles.matfile.geo.Nodes_ventricles;
    [m,im] = min(sum(abs(xyz - repmat(pos,[size(xyz,1) 1])),2));
    po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
    
    if handles.tag_menu.Value==1
        handles.matfile.Markers.dt(ib,im) = nan;
        po.FaceVertexCData = handles.matfile.Markers.dt(ib,:)'-handles.matfile.spikes(ib);
        
    elseif handles.tag_menu.Value==2
        handles.matfile.Markers.rt_Wyatt(ib,im) = nan;
        po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.spikes(ib);
        
    elseif handles.tag_menu.Value==3
        handles.matfile.Markers.dt(ib,im) = nan;
        handles.matfile.Markers.rt_Wyatt(ib,im) = nan;
        po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.Markers.dt(ib,im);
        
    elseif handles.tag_menu.Value>3
        return
    end
    
elseif isequal(get(info_struct.Target.Parent,'type'),'axes')
    pos = info_struct.Position;
    x = handles.matfile.Markers.rt_Wyatt(ib,:)-handles.matfile.spikes(1);
    ii = find(x==pos(1));
    jj = zeros(1,length(ii));
    for i = 1:length(ii)
        jj(i) = (handles.matfile.signals_proc(round(handles.matfile.Markers.rt_Wyatt(ib,ii(i))/1000*handles.matfile.ParamSig.frequency),ii(i))==pos(2));
    end
    ic = ii(logical(jj));
    handles.matfile.Markers.rt_Wyatt(ib,ic) = nan;
    if ~isempty(ic)
        po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
        po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.spikes(ib);
        handles.tag_menu.Value = 2;
    end
    clear ic;
    %
    x = handles.matfile.Markers.dt(ib,:)-handles.matfile.spikes(1);
    ii = find(x==pos(1));
    jj = zeros(1,length(ii));
    for i = 1:length(ii)
        jj(i) = (handles.matfile.signals_proc_AT(round(handles.matfile.Markers.dt(ib,ii(i))/1000*handles.matfile.ParamSig.frequency),ii(i))==pos(2));
    end
    handles.matfile.Markers.dt(ib,ii(find(jj))) = nan;
    
    ic = ii(logical(jj));
    if ~isempty(ic)
        po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
        po.FaceVertexCData = handles.matfile.Markers.dt(ib,:)'-handles.matfile.spikes(ib);
        handles.tag_menu.Value = 1;
    end
    clear ic;
    
    pp = findobj(handles.tag_egm,'marker','o','markersize',6);
    jj= sum([[pp.XData];[pp.YData]]==[pos(1);pos(2)]*ones(1,length(pp)))==2;
    delete(pp(jj));
    
    pp = findobj(handles.tag_egm,'marker','x','markersize',6);
    jj= sum([[pp.XData];[pp.YData]]==[pos(1);pos(2)]*ones(1,length(pp)))==2;
    delete(pp(jj));
    
    
end
guidata(hObject,handles)



% --- Executes on button press in tag_all_chan.
function tag_all_chan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_all_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_all_chan



function tag_plot_egm_numb_chan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_egm_numb_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_plot_egm_numb_chan as text
%        str2double(get(hObject,'String')) returns contents of tag_plot_egm_numb_chan as a double


% --- Executes during object creation, after setting all properties.
function tag_plot_egm_numb_chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_plot_egm_numb_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_markers.
function tag_markers_Callback(hObject, eventdata, handles)
% hObject    handle to tag_markers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Markers_input_ECGI(hObject)
handles = guidata(hObject);

fs = handles.matfile.ParamSig.frequency;
UEG = handles.matfile.signals;
DTmax = handles.MarkersInput.DTmax;
RTlimit = [handles.MarkersInput.min_RT handles.MarkersInput.max_RT];
RT_filter = handles.MarkersInput.RT_BW;
AT_filter = handles.MarkersInput.AT_BW;

[Markers] = ECGI_dt_rt_GUI(UEG,fs,DTmax,RTlimit,RT_filter,AT_filter);

handles.matfile.Markers = Markers;
guidata(hObject,handles)



function tag_beat_number_Callback(hObject, eventdata, handles)
% hObject    handle to tag_beat_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_beat_number as text
%        str2double(get(hObject,'String')) returns contents of tag_beat_number as a double


% --- Executes during object creation, after setting all properties.
function tag_beat_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_beat_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function tag_slider_beat_Callback(hObject, eventdata, handles)
% hObject    handle to tag_slider_beat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

delete(findobj(handles.tag_egm,'linewidth',2,'color','r','linestyle','--'));

ib = round(handles.tag_slider_beat.Value);
ib0 = str2double(handles.tag_beat_number.String);
if abs(ib-ib0)>1
    handles.tag_slider_beat.Value = ib0;
    return
end

handles.tag_beat_number.String = num2str(ib);
po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');

% if handles.tag_menu.Value==2
%     if ~isempty(po)
%         po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
%         cb = get(handles.tag_mesh,'clim');
%         po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.spikes(ib);
%         set(handles.tag_mesh,'clim',cb);
%     end
% end
% 
% if handles.tag_menu.Value==1
%     po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
%     if ~isempty(po)
%         cb = get(handles.tag_mesh,'clim');
%         
%         po.FaceVertexCData = handles.matfile.Markers.dt(ib,:)'-handles.matfile.spikes(ib);
%         set(handles.tag_mesh,'clim',cb);
%     end
% end

if ~isempty(po)
    cb = get(handles.tag_mesh,'clim');
    switch handles.tag_menu.Value
        case 2 % RT
            po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.spikes(ib);      
        case 1 % AT
            po.FaceVertexCData = handles.matfile.Markers.dt(ib,:)'-handles.matfile.spikes(ib);
        case 3 % ARI
            po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.Markers.dt(ib,:)';
    end
    set(handles.tag_mesh,'clim',cb);
end

if handles.tag_interpol.Value==1
    clear po
    po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
    
    v = po.FaceVertexCData;
    xyz = po.Vertices;
    iiok = ~isnan(v);
    F = scatteredInterpolant(xyz(iiok,1),xyz(iiok,2),xyz(iiok,3),v(iiok));
    v_int = F(xyz(:,1),xyz(:,2),xyz(:,3));
    po.FaceVertexCData = v_int;
    handles.map_values_ori = v;
    guidata(hObject,handles);
end
    
tag_plot_egm_Callback(hObject, eventdata, handles)




    
% pa = findobj(handles.tag_egm,'type','line');
% if ~isempty(pa)
%     if ib == length(handles.matfile.spikes)
%         xlim(handles.tag_egm,[handles.matfile.spikes(ib)-500 handles.matfile.spikes(ib)+1500]-handles.matfile.spikes(1));
%     else
%         xlim(handles.tag_egm,[handles.matfile.spikes(ib)-500 handles.matfile.spikes(ib+1)+500]-handles.matfile.spikes(1));
%     end
% end
% plot(handles.tag_egm,[handles.matfile.spikes(ib)]*[1 1]-handles.matfile.spikes(1),get(handles.tag_egm,'ylim'),'--r','linewidth',2)



% --- Executes during object creation, after setting all properties.
function tag_slider_beat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_slider_beat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in tag_Convex_Hull.
function tag_Convex_Hull_Callback(hObject, eventdata, handles)
% hObject    handle to tag_Convex_Hull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_Convex_Hull
if  get(handles.tag_Convex_Hull,'value')
    set(handles.tag_view_mesh,'value',0)
    set(handles.tag_view_mesh_interp,'value',0)
else
    set(handles.tag_view_mesh,'value',1)
end
% --- Executes on button press in tag_view_mesh.
function tag_view_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to tag_view_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_view_mesh

hp = findobj(handles.tag_mesh,'type','patch');
if get(handles.tag_view_mesh,'value')
    if ~isempty(hp)
        set(hp,'visible','on')
    end
    set(handles.tag_Convex_Hull,'value',0)
    set(handles.tag_view_mesh_interp,'value',0)
else
    if ~isempty(hp)
        set(hp,'visible','off')
    end
    set(handles.tag_Convex_Hull,'value',1)
end


% --- Executes on button press in tag_alpha.
function tag_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to tag_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_alpha
hp = findobj(handles.tag_mesh,'type','patch');
if ~isempty(hp)
    if get(handles.tag_alpha,'value')
        set(hp,'facealpha',0.25)
    else
        set(hp,'facealpha',1)
    end
end


% --- Executes on button press in tag_RVI.
function tag_RVI_Callback(hObject, eventdata, handles)
% hObject    handle to tag_RVI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ib = str2double(handles.tag_beat_number.String);
dt = handles.matfile.Markers.dt(ib,:) - handles.matfile.spikes(ib);
rt = handles.matfile.Markers.rt_Wyatt(ib,:) - handles.matfile.spikes(ib);
D = str2double(handles.tag_RVI_d.String);
rvi_type = handles.tag_RVI_type.Value;
xyz = handles.matfile.geo.Nodes_ventricles;
values = nan(length(dt),1);
Nnodes = nan(size(values));

if ~isfield(handles.matfile,'Disp_Local');
    handles.matfile.Disp_Local.range_AT = nan(length(handles.matfile.spikes),length(values));
    handles.matfile.Disp_Local.range_ARI = nan(length(handles.matfile.spikes),length(values));
    handles.matfile.Disp_Local.range_RT =nan(length(handles.matfile.spikes),length(values));
    handles.matfile.Disp_Local.range_AT_norm = nan(length(handles.matfile.spikes),length(values));
    handles.matfile.Disp_Local.range_ARI_norm = nan(length(handles.matfile.spikes),length(values));
    handles.matfile.Disp_Local.range_RT_norm = nan(length(handles.matfile.spikes),length(values));
end


for ic = 1:length(dt)
    d = sqrt(sum( (xyz  - repmat(xyz(ic,:),[size(xyz ,1),1])).^2 ,2));
    [Dim,ii] = sort(d,'ascend');
    ii(1)=[];Dim(1)=[];
    ii = ii(Dim<=D);
    d = d(ii);
    
    if ~isempty(ii)
        RVIall = rt(:,ic)*ones(1,length(ii)) - dt(:,ii);
        if rvi_type==1
            values(ic) =  nanmin(RVIall);
            Nnodes(ic) = sum(~isnan(RVIall));
        elseif rvi_type==2
            values(ic) =  nanmean(RVIall);
            Nnodes(ic) = sum(~isnan(RVIall));
        end
        
        %%
        handles.matfile.Disp_Local.range_AT(ib,ic) = max(abs(dt(ii)-dt(ic)));
        handles.matfile.Disp_Local.range_ARI(ib,ic) = max(abs(rt(ii)-rt(ic)));
        handles.matfile.Disp_Local.range_RT(ib,ic) = max(abs(rt(ii)-dt(ii)-(rt(ic)-dt(ic))));
        % -
        handles.matfile.Disp_Local.range_AT_norm(ib,ic) = max(abs(dt(ii)-dt(ic))./d.');
        handles.matfile.Disp_Local.range_ARI_norm(ib,ic) = max(abs(rt(ii)-rt(ic))./d.');
        handles.matfile.Disp_Local.range_RT_norm(ib,ic) = max(abs(rt(ii)-dt(ii)-(rt(ic)-dt(ic)))./d.');
        
    end
end
if ~isfield(handles.matfile,'RVI')
    handles.matfile.RVI = struct;
end

if rvi_type==1
    if ~isfield(handles.matfile.RVI,['min_D_',handles.tag_RVI_d.String]);
        handles.matfile.RVI.(['min_D_',handles.tag_RVI_d.String]) = nan(length(handles.matfile.spikes),size(handles.matfile.signals,2));
        handles.matfile.RVI.(['min_D_',handles.tag_RVI_d.String,'_Nnodes']) = nan(length(handles.matfile.spikes),size(handles.matfile.signals,2));
    end
    handles.matfile.RVI.(['min_D_',handles.tag_RVI_d.String])(ib,:) = values(:).';
    handles.matfile.RVI.(['min_D_',handles.tag_RVI_d.String,'_Nnodes'])(ib,:) = Nnodes(:).';
elseif rvi_type==2
    if ~isfield(handles.matfile.RVI,['avg_D_',handles.tag_RVI_d.String]);
        handles.matfile.RVI.(['avg_D_',handles.tag_RVI_d.String]) = nan(length(handles.matfile.spikes),size(handles.matfile.signals,2));
        handles.matfile.RVI.(['avg_D_',handles.tag_RVI_d.String,'_Nnodes']) = nan(length(handles.matfile.spikes),size(handles.matfile.signals,2));
    end
    handles.matfile.RVI.(['avg_D_',handles.tag_RVI_d.String])(ib,:) =values(:).';
    handles.matfile.RVI.(['avg_D_',handles.tag_RVI_d.String,'_Nnodes'])(ib,:) = Nnodes(:).';
end

po = findobj(handles.figure1,'type','patch');
po = findobj(po,'facecolor','interp');


values_int  = values;
set(po,'FaceVertexCData',values_int)
set(get(po,'parent'),'clim',[prctile(values_int,5) prctile(values_int,95)])
pc = findobj(handles.figure1,'type','colorbar');
set(pc,'ylim',[prctile(values_int,5) prctile(values_int,95)])
if rvi_type==1
    title(pc,['RVI_{min} [D=',num2str(D),'mm]'])
elseif rvi_type==2
    title(pc,['RVI_{avg} [D=',num2str(D),'mm]'])
end
handles.tag_menu.Value = 6;
guidata(hObject,handles);

% --- Executes on selection change in tag_RVI_type.
function tag_RVI_type_Callback(hObject, eventdata, handles)
% hObject    handle to tag_RVI_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_RVI_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_RVI_type

% --- Executes during object creation, after setting all properties.
function tag_RVI_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_RVI_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_RVI_d_Callback(hObject, eventdata, handles)
% hObject    handle to tag_RVI_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_RVI_d as text
%        str2double(get(hObject,'String')) returns contents of tag_RVI_d as a double


% --- Executes during object creation, after setting all properties.
function tag_RVI_d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_RVI_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_view_mesh_interp.
function tag_view_mesh_interp_Callback(hObject, eventdata, handles)
% hObject    handle to tag_view_mesh_interp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_view_mesh_interp
if  get(handles.tag_view_mesh_interp,'value')
    set(handles.tag_view_mesh,'value',0)
    set(handles.tag_Convex_Hull,'value',0)
else
    set(handles.tag_view_mesh_interp,'value',1)
end


% --- Executes on button press in tag_modify_ALL.
function tag_modify_ALL_Callback(hObject, eventdata, handles)
% hObject    handle to tag_modify_ALL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% GUI_EGM_map_viewer_Ensite_modify_all_info(hObject)

Output = GUI_EGM_map_viewer_Ensite_modify_all_info(hObject);

%%
if ~isempty(Output)
    mode_type = Output.do_type;
    do_AT = Output.do_AT;
    do_RT = Output.do_RT;
    
    if do_AT==0&do_RT==0
        msgbox('Select AT or RT or both');
        return
    end
    AT_min = Output.tab{1,1};
    AT_max = Output.tab{1,2};
    RT_min = Output.tab{2,1};
    RT_max = Output.tab{2,2};
    Fcut_AT = Output.tab{1,3};
    Fcut_RT = Output.tab{2,3};
    
    
    % ib = str2double(handles.tag_beat_number.String);
    ib_tot = Output.beats;
    if length(ib_tot)>1
        hb=waitbar(0/length(ib_tot),'Finding markers ...');
    end
    for iter = 1:length(ib_tot);
        
        ib = ib_tot(iter);
        if do_AT
            if mode_type==1
                ic = 1:size(handles.matfile.signals,2);
            else
                ic = find(isnan(handles.matfile.Markers.dt(ib,:)));
            end
            
            %             W = round(AT_min/1000*handles.matfile.ParamSig.frequency)+1 : round([AT_max]/1000*handles.matfile.ParamSig.frequency);
            %             W(W<1)=1;
            %             W(W>size(handles.matfile.signals,1))=size(handles.matfile.signals,1);
            
            %             if ib < length(handles.matfile.spikes);
            %                  Ws = round(handles.matfile.spikes(ib)/1000*handles.matfile.ParamSig.frequency):round(handles.matfile.spikes(ib+1)/1000*handles.matfile.ParamSig.frequency);
            %             else
            %                  Ws = round(handles.matfile.spikes(ib)/1000*handles.matfile.ParamSig.frequency):size(handles.matfile.signals,1);
            %             end
            
            Ws = round( (handles.matfile.spikes(ib)+AT_min)/1000*handles.matfile.ParamSig.frequency) : round( (handles.matfile.spikes(ib)+AT_max)/1000*handles.matfile.ParamSig.frequency);
            
            Ws(Ws<1)=[];
            %             t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);
            
            [sig2,~] = replace_spikes_mo(handles.matfile.signals,round(handles.matfile.spikes/1000*1200),25,10);
            f0 = Fcut_AT; %Hz
            [b, a] = butter(3,f0/(handles.matfile.ParamSig.frequency/2), 'low');
            signals_proc_AT = filtfilt(b,a,sig2);
            handles.matfile.signals_proc_Ensite_AT = signals_proc_AT;
            X = signals_proc_AT(Ws,ic);
            
            
            Xd = diff(X);Xd2 = diff(Xd);Xd3 = diff(Xd2);
            ii = Xd2([1:end-1],:).*Xd2([2 : end],:)<0 & Xd([1:end-2],:)<0 & Xd3([1: end],:)>0; % zeros of II derivative & I derivative negative & III derivative > 0
            for j = 1:size(ii,2)
                tt = find(ii(:,j));
                if ~isempty(tt)
                    [~,kk]=min(Xd(tt,j));
                    handles.matfile.Markers.dt(ib,ic(j)) = [tt(kk)+Ws(1)]/handles.matfile.ParamSig.frequency*1000;
                end
            end
            
            
            if handles.tag_menu.Value==1
                po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
                if ~isempty(po)
                    cb = get(handles.tag_mesh,'clim');
                    
                    po.FaceVertexCData = handles.matfile.Markers.dt(ib,:)'-handles.matfile.spikes(ib);
                    set(handles.tag_mesh,'clim',cb);
                end
            end
        end
        
        if do_RT
            
            if mode_type==1
                ic = 1:size(handles.matfile.signals,2);
            else
                ic = find(isnan(handles.matfile.Markers.rt_Wyatt(ib,:)));
            end
            
            
            if ib < length(handles.matfile.spikes);
                Ws = round(handles.matfile.spikes(ib)/1000*handles.matfile.ParamSig.frequency):round(handles.matfile.spikes(ib+1)/1000*handles.matfile.ParamSig.frequency);
            else
                Ws = round(handles.matfile.spikes(ib)/1000*handles.matfile.ParamSig.frequency):size(handles.matfile.signals,1);
            end
            Ws(Ws<1)=[];
            t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);
            X = handles.matfile.signals(Ws,ic);
            
            W = round(RT_min/1000*handles.matfile.ParamSig.frequency)+1 : round([RT_max]/1000*handles.matfile.ParamSig.frequency);
            W(W<1)=1;
            W(W>size(X,1)-2) = [];
            if isempty(W)
                continue
            end
            
            f0 = Fcut_RT; %Hz
            [b, a] = butter(3,f0/(handles.matfile.ParamSig.frequency/2), 'low');
            signals_proc_RT = filtfilt(b,a,handles.matfile.signals);
            handles.matfile.signals_proc_Ensite_RT = signals_proc_RT;
            X = signals_proc_RT(Ws,ic);
            
            Xd = diff(X);Xd2 = diff(Xd);Xd3 = diff(Xd2);
            
            ii = Xd2([W(1):W(end)-1],:).*Xd2([W(1)+1 : W(end)],:)<0 & Xd([W(1): W(end)-1],:)>0 & Xd3([W(1): W(end)-1],:)<=0;
            for j = 1:size(ii,2)
                tt = find(ii(:,j));
                
                if ~isempty(tt)
                    [~,kk]=max(Xd(tt,j));
                    handles.matfile.Markers.rt_Wyatt(ib,ic(j)) = [tt(kk)+W(1)+Ws(1)-1]/handles.matfile.ParamSig.frequency*1000;
                    
                end
            end
            
            %
            if handles.tag_menu.Value==2
                po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
                if ~isempty(po)
                    po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
                    cb = get(handles.tag_mesh,'clim');
                    po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.spikes(ib);
                    set(handles.tag_mesh,'clim',cb);
                end
            end
        end
        if length(ib_tot)>1
            waitbar(iter/length(ib_tot),hb);
        end
    end
    
    %Plot beat of main window
    if length(ib_tot)>1
        ib = str2double(handles.tag_beat_number.String);
        if handles.tag_menu.Value==2
            if ~isempty(po)
                po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
                cb = get(handles.tag_mesh,'clim');
                po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.spikes(ib);
                set(handles.tag_mesh,'clim',cb);
            end
        end
        
        if handles.tag_menu.Value==1
            po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
            if ~isempty(po)
                cb = get(handles.tag_mesh,'clim');
                
                po.FaceVertexCData = handles.matfile.Markers.dt(ib,:)'-handles.matfile.spikes(ib);
                set(handles.tag_mesh,'clim',cb);
            end
        end
    end
    
    guidata(hObject,handles)
    delete(findall(gcf,'Type','hggroup'));
    if length(ib_tot)>1
        delete(hb);
    end
end


% --- Executes on button press in tag_done.
function tag_done_Callback(hObject, eventdata, handles)
% hObject    handle to tag_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles_main = guidata(handles.hObject_main);
handles_main.MarkersC = handles.matfile.Markers;
if isfield(handles.matfile,'RVI')
    handles_main.RVI = handles.matfile.RVI;
end
if isfield(handles.matfile,'signals_proc_Ensite_RT')
    handles_main.signals_proc_Ensite_RT = handles.matfile.signals_proc_Ensite_RT;
end
if isfield(handles.matfile,'signals_proc_Ensite_AT')
    handles_main.signals_proc_Ensite_AT = handles.matfile.signals_proc_Ensite_AT;
end

if isfield(handles.matfile,'geo')
    handles_main.geo = handles.matfile.geo;
end
guidata(handles.hObject_main,handles_main);
close(handles.figure1)


% --- Executes on button press in tag_fig_summary.
function tag_fig_summary_Callback(hObject, eventdata, handles)
% hObject    handle to tag_fig_summary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


PROMPT = 'Enter index of selected VT beat';
NAME = 'Info';
NUMLINES = 1;
DEFAULTANSWER = {num2str(length(handles.matfile.spikes) - 1)};
ANSWER = inputdlg(PROMPT,NAME,NUMLINES,DEFAULTANSWER);

ibVT = str2double(ANSWER{1});
ib = str2double(handles.tag_beat_number.String);
if handles.tag_RVI_type.Value == 1
    name = ['min_D_',handles.tag_RVI_d.String];
else
    name = ['avg_D_',handles.tag_RVI_d.String];
end

if ~isfield(handles.matfile,'RVI');
    msgbox(['Calculate RVI-',name,' for beat=',num2str(ib)]);
    return
end

fig_maps = figure;
set(fig_maps,'units','centimeters','paperunits','centimeters','position',[3 3 20 8],...,
    'paperposition',[0 0 20 8],'papersize',[20 8])

%% Mesh

connectivity = handles.matfile.geo.Mesh_ventricles;
values = nan(size(handles.matfile.Markers.dt,2),4);
values(:,1) = handles.matfile.Markers.dt(ib,:) - handles.matfile.spikes(ib);
values(:,2) = handles.matfile.Markers.rt_Wyatt(ib,:) - handles.matfile.spikes(ib);
values(:,3) = handles.matfile.RVI.(name)(ib,:);
values(:,4) = handles.matfile.Markers.dt(ibVT,:) - handles.matfile.spikes(ibVT);

colorLabel{1} = 'AT (ms)';
colorLabel{2} = 'RT (ms)';
if handles.tag_RVI_type.Value == 1
    colorLabel{3} = ['RVI_min [D=',handles.tag_RVI_d.String,'] (ms)'];
else
    colorLabel{3} = ['RVI_avg [D=',handles.tag_RVI_type.String,'] (ms)'];
end
colorLabel{4} = 'AT(VT beat) (ms)';

for i = 1:4
    MV.faces = connectivity;
    MV.vertices =  handles.matfile.geo.Nodes_ventricles;
    MV.faceVertexCData = values(:,i);
    
    ax(i) = subplot(1,4,i);
    meshObject(i) = patch(MV,'parent',ax(i));
    cbar(i) = colorbar;
    set(ax(i),'clim',[prctile(values(:,i),2) prctile(values(:,i),98)])
    set(cbar(i),'Location','southoutside')
    title(cbar(i),colorLabel{i},'fontsize',12);
end




%     set(cbar,'Location','southoutside')
%     title(cbar,colorLabel,'fontsize',10);
%     set(cbar,'position',[.02 .05 .4 .025],'xlim',Clim)
%     % camorbit(50,-10);
rotate3d(handles.tag_mesh);

if exist([pwd,'\GUI_egm_mFiles\Geo_Chann\colormap_carto.mat'],'file')>0
    load([pwd,'\GUI_egm_mFiles\Geo_Chann\colormap_carto.mat']);
    cmap = cmap([1 2:2:end],:);
    colormap(cmap)
else
    cmap = jet;
    cmap = cmap(end:-1:1,:);
    colormap(cmap)
end

% set(meshObject,'markersize',2,'marker','o','markerFaceColor',[1 1 1]*.5)

if isfield(handles.matfile.geo,'Mesh_Valves')
    MV_valves.faces = handles.matfile.geo.Mesh_Valves;
    MV_valves.vertices =  handles.matfile.geo.Nodes_Valves;
    MV_valves.faceVertexCData = 1;
else
    MV_valves = [];
end
%
if isfield(handles.matfile.geo,'Mesh_LAD')
    MV_lad.faces = handles.matfile.geo.Mesh_LAD;
    MV_lad.vertices =  handles.matfile.geo.Nodes_LAD;
    MV_lad.faceVertexCData = 1;
else
    MV_lad = [];
end
%

% LAD
if isstruct(MV_lad)
    
    for i = 1:4
        hold(ax(i),'on')
        meshObject_LAD(i) = patch(MV_lad,'parent',ax(i));
        set(meshObject_LAD(i),'edgecolor','none','facecolor',[1 1 1]*0.3);
    end
end

% Valves
if isstruct(MV_valves)
    for i = 1:4
        hold(ax(i),'on')
        meshObject_valves(i) = patch(MV_valves,'parent',ax(i));
        set(meshObject_valves(i),'edgecolor','none','facecolor',[1 1 1]*0.3);
    end
end

L = .22;DL = 0.01;L0=0.02;
H=.85;
set(ax(1),'position',[L0+(L+DL)*0 .15 L H])
set(ax(2),'position',[L0+(L+DL)*1 .15 L H])
set(ax(3),'position',[L0+(L+DL)*2 .15 L H])
set(ax(4),'position',[L0+(L+DL)*3+DL*4 .15 L H])

set(cbar(1),'position',[L0+(L+DL)*0 .10 L .025])
set(cbar(2),'position',[L0+(L+DL)*1 .1 L .025])
set(cbar(3),'position',[L0+(L+DL)*2 .1 L .025])
set(cbar(4),'position',[L0+(L+DL)*3+DL*4 .1 L .025])

xyz = handles.matfile.geo.Nodes_ventricles;
set(ax,'xlim',[min(xyz(:,1)) max(xyz(:,1))],'ylim',[min(xyz(:,2)) max(xyz(:,2))],'zlim',[min(xyz(:,3)) max(xyz(:,3))])
axis(ax,'off','image')
h_link = linkprop(ax,'view');
setappdata(fig_maps,'x_axis_linkprop',h_link);
set(meshObject,'edgecolor','none','facecolor','interp','facelighting','gouraud',...
    'ambientstrength',0.9,'specularcolorreflectance',0.4,'specularexponent',3,'specularstrength',0.3);
set(ax,'fontsize',11)
a1 = annotation('textbox',[L0 .9 L0+(L+DL)*2+L .08],'string','RVI BEAT');
a2 = annotation('textbox',[L0+(L+DL)*3+DL*4 .9 L .08],'string','VT BEAT');
set([a1 a2],'horizontalalignment','center','fontsize',13,'verticalalignment','middle')


% --- Executes on selection change in tag_find_chan.
function tag_find_chan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_find_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_find_chan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_find_chan

cc = findobj(handles.tag_mesh,'marker','x');
if ~isempty(cc)&~handles.tag_hold_on.Value
    delete(cc)
end
delete(findall(gcf,'Type','hggroup'))
datacursormode off

ic = handles.tag_find_chan.Value;
hold(handles.tag_mesh,'on');
aa = plot3(handles.matfile.geo.Nodes_ventricles(ic,1),handles.matfile.geo.Nodes_ventricles(ic,2),handles.matfile.geo.Nodes_ventricles(ic,3),'x');
set(aa,'markerfacecolor',[1 1 1]*.6,'markeredgecolor','w','markersize',20,'linewidth',2)


% --- Executes during object creation, after setting all properties.
function tag_find_chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_find_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_edit_CBmin_Callback(hObject, eventdata, handles)
% hObject    handle to tag_edit_CBmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_edit_CBmin as text
%        str2double(get(hObject,'String')) returns contents of tag_edit_CBmin as a double


% --- Executes during object creation, after setting all properties.
function tag_edit_CBmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_edit_CBmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tag_edit_CBmax_Callback(hObject, eventdata, handles)
% hObject    handle to tag_edit_CBmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_edit_CBmax as text
%        str2double(get(hObject,'String')) returns contents of tag_edit_CBmax as a double


% --- Executes during object creation, after setting all properties.
function tag_edit_CBmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_edit_CBmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_button_CBreset.
function tag_button_CBreset_Callback(hObject, eventdata, handles)
% hObject    handle to tag_button_CBreset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Clim(1) = str2double(handles.tag_edit_CBmin.String);
Clim(2) = str2double(handles.tag_edit_CBmax.String);

pcbar = findobj(handles.tag_mesh,'tag','colorbar');
set(handles.tag_mesh,'Clim',Clim);
set(handles.tag_mesh,'Clim',Clim);
cbar = findobj(handles.figure1,'tag','Colorbar');
set(cbar,'xlim',Clim)


% --- Executes during object creation, after setting all properties.
function tag_r_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_r_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in tag_tab_print.
function tag_tab_print_Callback(hObject, eventdata, handles)
% hObject    handle to tag_tab_print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fname = handles.matfile.ParamSig.name_file;
ii = find(fname == '.');
fname(ii(end):end) = [];
[FILENAME, PATHNAME] = uiputfile([fname,'.xls'],'Select folder where to save the spreadsheet');
tabname = [PATHNAME,FILENAME];

xyz = handles.matfile.ParamSig.Virtuals_geo.xyz;
tab_geo = [[1:size(xyz,1)]',xyz];
Dtot = [10 15 20];
HL = 2*length(handles.matfile.spikes)+1;
H = waitbar(1/HL,'Printing ...');

% -
Disp_Local.range_AT = nan(length(handles.matfile.spikes),size(xyz,1),length(Dtot));
Disp_Local.range_ARI = nan(length(handles.matfile.spikes),size(xyz,1),length(Dtot));
Disp_Local.range_RT =nan(length(handles.matfile.spikes),size(xyz,1),length(Dtot));
Disp_Local.range_AT_norm = nan(length(handles.matfile.spikes),size(xyz,1),length(Dtot));
Disp_Local.range_ARI_norm = nan(length(handles.matfile.spikes),size(xyz,1),length(Dtot));
Disp_Local.range_RT_norm = nan(length(handles.matfile.spikes),size(xyz,1),length(Dtot));
RVImin = nan(length(handles.matfile.spikes),size(xyz,1),length(Dtot));
RVImean = nan(length(handles.matfile.spikes),size(xyz,1),length(Dtot));
dt_all = nan(length(handles.matfile.spikes),size(xyz,1));
rt_all = nan(length(handles.matfile.spikes),size(xyz,1));

% -
for ib = 1 : length(handles.matfile.spikes)
    waitbar((1+ib)/HL,H);
    dt = handles.matfile.Markers.dt(ib,:) - handles.matfile.spikes(ib);
    rt = handles.matfile.Markers.rt_Wyatt(ib,:) - handles.matfile.spikes(ib);
    dt_all(ib,:) = dt;
    rt_all(ib,:) = rt;
    
    for id = 1:length(Dtot)
        D = Dtot(id);
        
        for ic = 1:length(dt)
            d = sqrt(sum( (xyz  - repmat(xyz(ic,:),[size(xyz ,1),1])).^2 ,2));
            [Dim,ii] = sort(d,'ascend');
            ii(1)=[];Dim(1)=[];
            ii = ii(Dim<=D);
            d = d(ii);
            
            if ~isempty(ii)
                RVIall = rt(:,ic)*ones(1,length(ii)) - dt(:,ii);
                %%
                Disp_Local.range_AT(ib,ic,id) = max(abs(dt(ii)-dt(ic)));
                Disp_Local.range_ARI(ib,ic,id) = max(abs(rt(ii)-rt(ic)));
                Disp_Local.range_RT(ib,ic,id) = max(abs(rt(ii)-dt(ii)-(rt(ic)-dt(ic))));
                % -
                Disp_Local.range_AT_norm(ib,ic,id) = max(abs(dt(ii)-dt(ic))./d.');
                Disp_Local.range_ARI_norm(ib,ic,id) = max(abs(rt(ii)-rt(ic))./d.');
                Disp_Local.range_RT_norm(ib,ic,id) = max(abs(rt(ii)-dt(ii)-(rt(ic)-dt(ic)))./d.');
                % -
                RVImin(ib,ic,id) = min(RVIall);
                RVImean(ib,ic,id) = nanmean(RVIall);
            end
        end
        
        
    end
    
end

tab_RVI = {'RVImin','RVImean','RVIminC','RVImeanC','DAT','DARI','DRT','DAT(N)','DARI(N)','DRT(N)',' '};
tab_header = [{'chan','x','y','z',' ','AT','ARI','RT',' '} repmat(tab_RVI,1,length(Dtot))];

for ib = 1 : length(handles.matfile.spikes)
    
    tab_markers = [dt_all(ib,:)' rt_all(ib,:)'-dt_all(ib,:)' rt_all(ib,:)'];
    
    tab_RVI = [];tab_header_sup = [];
    for id = 1:length(Dtot)
        %         tab = [RVImin(ib,:,id)' RVImean(ib,:,id)' Disp_Local.range_AT(ib,:,id)' Disp_Local.range_ARI(ib,:,id)' Disp_Local.range_RT(ib,:,id)'   Disp_Local.range_AT_norm(ib,:,id)' Disp_Local.range_ARI_norm(ib,:,id)' Disp_Local.range_RT_norm(ib,:,id)' nan(size(xyz,1),1)];
        tab = [RVImin(ib,:,id)' RVImean(ib,:,id)' RVImin(ib,:,id)'-nanmedian(RVImin(ib,:,id)') RVImean(ib,:,id)'-nanmedian(RVImean(ib,:,id)') Disp_Local.range_AT(ib,:,id)' Disp_Local.range_ARI(ib,:,id)' Disp_Local.range_RT(ib,:,id)'   Disp_Local.range_AT_norm(ib,:,id)' Disp_Local.range_ARI_norm(ib,:,id)' Disp_Local.range_RT_norm(ib,:,id)' nan(size(xyz,1),1)];
        
        tab_RVI = [tab_RVI tab];
        tab_header_sup = [tab_header_sup [{['D=',num2str(Dtot(id)),' mm']},cell(1,size(tab,2)-1)]];
    end
    
    tab_beat = [tab_geo nan(size(xyz,1),1) tab_markers nan(size(xyz,1),1) tab_RVI];
    tab_beat_cell = [[{'','Geometry','','','','Markers','','',''},tab_header_sup];tab_header;num2cell(tab_beat)];
    
    
    waitbar((1+length(handles.matfile.spikes)+ib)/HL,H);
    sheetname = ['b=',num2str(ib)];
    xlswrite(tabname,tab_beat_cell,sheetname,'A1');
    
    %%
    tab_beat_stat.dist = prctile(tab_beat,[5 10 25 50 75 90 95]);
    tab_beat_stat.values(1,:) = nanmean(tab_beat);
    tab_beat_stat.values(2,:) = nanstd(tab_beat);
    tab_beat_stat.values(3,:) = nanmedian(tab_beat);
    tab_beat_stat.values(4,:) = mad(tab_beat,1);
    tab_beat_stat.values(5,:) = range(tab_beat);
    tab_beat_stat.values(6,:) = prctile(tab_beat,95)-prctile(tab_beat,5);
    
    tab_beat_stat_header1 = {'5th','10th','25th','50th','75th','90th','95th'};
    tab_beat_stat_header2 = {'mean','SD','median','mad','range','5-95th'};
    
    xlswrite(tabname,num2cell(tab_beat_stat.dist),sheetname,'A270');
    xlswrite(tabname,tab_beat_stat_header1(:),sheetname,'A270');
    
    xlswrite(tabname,num2cell(tab_beat_stat.values),sheetname,'A280');
    xlswrite(tabname,tab_beat_stat_header2(:),sheetname,'A280');
    
end

close(H)


% --- Executes on selection change in tag_export_fig.
function tag_export_fig_Callback(hObject, eventdata, handles)
% hObject    handle to tag_export_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_export_fig contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_export_fig


if handles.tag_export_fig.Value==5 % RVI+VT
    tag_fig_summary_Callback(hObject, eventdata, handles)
    handles.tag_export_fig.Value=1;
    
elseif handles.tag_export_fig.Value==2 % Just export the current map
    h_fig = figure;
    set(gcf,'units','centimeters','paperUnits','centimeters','position',[2 2 14 12],'paperposition',[0 0 14 12],'papersize',[14 12])
    
    new_handle = copyobj(handles.tag_mesh,h_fig,'legacy');
    axis(new_handle,'image')
    set(new_handle,'position',[0.02 0.05 0.85 0.9])
    cmap = colormap(handles.tag_mesh);
    colormap(new_handle,cmap);
    cbar = colorbar(new_handle);
    hcbar = findobj(handles.figure1,'type','colorbar');
    bcbar_title = get(hcbar(end),'title');
    set(cbar,'position',[0.87 0.05 0.04 0.8])
    
    title(cbar,bcbar_title.String)
    set(new_handle,'fontsize',12)
    handles.tag_export_fig.Value=1;
    
elseif handles.tag_export_fig.Value==3 % Just export the current map (x2)
    h_fig = figure;
    set(gcf,'units','centimeters','paperUnits','centimeters','position',[2 2 18 10],'paperposition',[0 0 18 10],'papersize',[18 10])
    
    new_handle(1) = copyobj(handles.tag_mesh,h_fig,'legacy');
    new_handle(2) = copyobj(handles.tag_mesh,h_fig,'legacy');
    
    axis(new_handle,'image')
    set(new_handle(1),'position',[0.02 0.05 0.45 0.9]);
    set(new_handle(2),'position',[0.47 0.05 0.45 0.9]);
    
    cmap = colormap(handles.tag_mesh);
    colormap(new_handle(1),cmap);
    colormap(new_handle(2),cmap);
    cbar = colorbar(new_handle(2));
    hcbar = findobj(handles.figure1,'type','colorbar');
    bcbar_title = get(hcbar(end),'title');
    set(cbar,'position',[0.9 0.1 0.02 0.8])
    
    title(cbar,bcbar_title.String)
    set(new_handle,'fontsize',12)
    handles.tag_export_fig.Value=1;
    
    
elseif handles.tag_export_fig.Value==4 % Just export the current map (x4)
    h_fig = figure;
    set(gcf,'units','centimeters','paperUnits','centimeters','position',[2 2 18 7],'paperposition',[0 0 18 7],'papersize',[18 7])
    
    new_handle(1) = copyobj(handles.tag_mesh,h_fig,'legacy');
    new_handle(2) = copyobj(handles.tag_mesh,h_fig,'legacy');
    new_handle(3) = copyobj(handles.tag_mesh,h_fig,'legacy');
    
    axis(new_handle,'image')
    
    set(new_handle(1),'position',[0.02 0.05 0.30 0.9]);
    set(new_handle(2),'position',[0.32 0.05 0.30 0.9]);
    set(new_handle(3),'position',[0.62 0.05 0.30 0.9]);
    
    cmap = colormap(handles.tag_mesh);
    colormap(new_handle(1),cmap);
    colormap(new_handle(2),cmap);
    colormap(new_handle(3),cmap);
    
    cbar = colorbar(new_handle(2));
    hcbar = findobj(handles.figure1,'type','colorbar');
    bcbar_title = get(hcbar(end),'title');
    set(cbar,'position',[0.9 0.1 0.02 0.8])
    
    title(cbar,bcbar_title.String)
    set(new_handle,'fontsize',12)
    handles.tag_export_fig.Value=1;
    
    
elseif handles.tag_export_fig.Value==4
    h_fig = figure;
    set(gcf,'units','centimeters','paperUnits','centimeters','position',[2 2 14 12],'paperposition',[0 0 14 12],'papersize',[14 12])
    
    
    hleg = findobj(handles.figure1,'type','legend');
    new_handle = copyobj([handles.tag_egm,hleg],h_fig,'legacy');
    set(new_handle(1),'position',[0.10 .12 0.85 0.8])
    
    set(new_handle,'fontsize',12)
    handles.tag_export_fig.Value=1;
end


% --- Executes during object creation, after setting all properties.
function tag_export_fig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_export_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_interpol.
function tag_interpol_Callback(hObject, eventdata, handles)
% hObject    handle to tag_interpol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_interpol

if handles.tag_interpol.Value;
    po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
    v = po.FaceVertexCData;
    xyz = po.Vertices;
    iiok = ~isnan(v);
    F = scatteredInterpolant(xyz(iiok,1),xyz(iiok,2),xyz(iiok,3),v(iiok));
    v_int = F(xyz(:,1),xyz(:,2),xyz(:,3));
    po.FaceVertexCData = v_int;
    handles.map_values_ori = v;
    guidata(hObject,handles);
else
    po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
    po.FaceVertexCData = handles.map_values_ori;
end


% --- Executes on button press in tag_plot_filt_RT.
function tag_plot_filt_RT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_filt_RT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_filt_RT
handles.tag_plot_filt_AT.Value = 0;
tag_plot_egm_Callback(hObject, eventdata, handles);
handles.tag_modify_type.Value=2;


% --- Executes on button press in tag_segment.
function tag_segment_Callback(hObject, eventdata, handles)
% hObject    handle to tag_segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_segment_mesh(handles);


% --- Executes on button press in tag_plot_filt_AT.
function tag_plot_filt_AT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_filt_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_filt_AT
handles.tag_plot_filt_RT.Value = 0;
tag_plot_egm_Callback(hObject, eventdata, handles);
handles.tag_modify_type.Value=1;


% --- Executes on button press in tag_include_NaNs_plot.
function tag_include_NaNs_plot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_include_NaNs_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_include_NaNs_plot
tag_plot_egm_Callback(hObject, [], handles)
