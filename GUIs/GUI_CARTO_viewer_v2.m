function varargout = GUI_CARTO_viewer_v2(varargin)
% GUI_CARTO_VIEWER_v2 MATLAB code for GUI_CARTO_viewer_v2.fig
%      GUI_CARTO_VIEWER_v2, by itself, creates a new GUI_CARTO_VIEWER_v2 or raises the existing
%      singleton*.
%
%      H = GUI_CARTO_VIEWER_v2 returns the handle to a new GUI_CARTO_VIEWER_v2 or the handle to
%      the existing singleton*.
%
%      GUI_CARTO_VIEWER_v2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CARTO_VIEWER_v2.M with the given input arguments.
%
%      GUI_CARTO_VIEWER_v2('Propguideerty','Value',...) creates a new GUI_CARTO_VIEWER_v2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_CARTO_viewer_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_CARTO_viewer_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_CARTO_viewer_v2

% Last Modified by GUIDE v2.5 28-Aug-2018 16:13:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_CARTO_viewer_v2_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_CARTO_viewer_v2_OutputFcn, ...
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


% --- Executes just before GUI_CARTO_viewer_v2 is made visible.
function GUI_CARTO_viewer_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_CARTO_viewer_v2 (see VARARGIN)

% Choose default command line output for GUI_CARTO_viewer_v2
handles.output = hObject;

%%
h = guidata(varargin{1});
handles.matfile.signals_proc = h.signals_proc;
handles.matfile.signals = h.signals;
if isfield(h,'signals_proc_AT')
    handles.matfile.signals_proc_AT = h.signals_proc_AT;
else
    handles.tag_plot_filt_AT.Visible = 0;
end
handles.matfile.spikes = h.spikes;

if ~isfield(h,'SNR');
    P = abs(fft(detrend(handles.matfile.signals,'constant'))).^2;
    ff= [1:size(P,1)]/size(P,1)*h.ParamSig.frequency;
    handles.matfile.SNR = 10*log10( sum(P(ff>2&ff<=40,:))./sum(P(ff>40&ff<=100,:)));
    clear P ff
else
    handles.matfile.SNR = h.SNR;
end

if isfield(h,'MarkersC')
    handles.matfile.Markers = h.MarkersC;
else
    handles.matfile.Markers = h.Markers;
end
handles.matfile.ParamSig = h.ParamSig;
%
if isfield(h,'RVI')
    handles.matfile.RVI = h.RVI;
end
%
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

% UIWAIT makes GUI_CARTO_viewer_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_CARTO_viewer_v2_OutputFcn(hObject, eventdata, handles)
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
function tag_plot_egm_Callback(hObject,i_chan, handles)
% hObject    handle to tag_plot_egm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xyz = handles.matfile.geo.xyz;

% these are channels that have delated and have all nans
cc = findobj(handles.tag_mesh,'marker','o','markersize',8);
cc = cc(end:-1:1);
iiko = find(strcmp({cc.Visible},'off'));
xyz(iiko,:) = nan;
clear iiko

cc = findobj(handles.tag_egm,'linestyle','--','color','k');
delete(cc);
dcm_obj = datacursormode(handles.figure1);
info_struct = getCursorInfo(dcm_obj);

if ~isnumeric(i_chan)
    if ~isempty(info_struct)
        pos = info_struct.Position;
    else
        ic = handles.tag_find_chan.Value;
        pos = xyz(ic,:);
    end
else
    pos = xyz(i_chan,:); % from delate
end

d = sqrt(sum( (xyz  - repmat(pos,[size(xyz ,1),1])).^2 ,2));
[~,im] = sort(d,'ascend');
im = im(1:str2double(get(handles.tag_plot_egm_numb_chan,'string')));


set(handles.tag_find_chan,'value',im(1));
handles.tag_r_plot.String = num2str(max(d(im)),'%1.1f');

t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);

%
if handles.tag_plot_filt.Value
    S = handles.matfile.signals_proc(:,im);
elseif handles.tag_plot_filt_raw.Value
    S = handles.matfile.signals(:,im);
elseif handles.tag_plot_filt_AT.Value
    S = handles.matfile.signals_proc_AT(:,im);
end

if handles.tag_median_map.Value
    ib = [1:length(handles.matfile.spikes)];
else
    ib = str2double(get(handles.tag_beat_number,'string'));
end

yy_lim = [];
xx_lim = [];


if get(handles.tag_hold_on,'value')==1;
    hold(handles.tag_egm,'on'),
    bb = findobj(handles.tag_egm,'tag','UEG');
    xx_lim = get(handles.tag_egm,'xlim')
else
    if isnumeric(i_chan)
        xx_lim = get(handles.tag_egm,'xlim');
        yy_lim = get(handles.tag_egm,'ylim');
    else
        xx_lim = [handles.matfile.spikes(ib)-100 handles.matfile.spikes(ib)+1000]-handles.matfile.spikes(1);
    end
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

strleg = cell(size(aueg));
for i= 1:length(aueg)
    strleg(i) = {['node=',num2str(im(i))]};
end
legend([aueg;bb],[strleg(:);bbl(:)])
set(aueg,'tag','UEG')

% legend(aueg,strleg)
hold(handles.tag_egm,'on'),
set(handles.tag_egm,'box','off')
xlabel(handles.tag_egm,'Time (ms)')

dt = handles.matfile.Markers.dt(ib,im) - handles.matfile.spikes(1);
rt = handles.matfile.Markers.rt_Wyatt(ib,im) - handles.matfile.spikes(1);
for j = 1:length(im)
    % AT
    %     [~,dtp] = find(t==dt(:,j));
    [~,dtp] = find([ones(length(dt(:,j)),1)*t - dt(:,j)*ones(1,length(t))]==0);
    plot(handles.tag_egm,t(dtp)',S(dtp,j),'marker','o','markerfacecolor',[1 1 1]*.6,'markeredgecolor','k','linestyle','none')
    % RT
    %     [~,rtp] = find(t==rt(:,j));
    [~,rtp] = find([ones(length(dt(:,j)),1)*t - rt(:,j)*ones(1,length(t))]==0);
    plot(handles.tag_egm,t(rtp)',S(rtp,j),'marker','x','markerfacecolor',[1 1 1]*.8,'markeredgecolor','k','linestyle','none')
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

if numel(ib)==1
    plot(handles.tag_egm,[handles.matfile.spikes(ib)]*[1 1]-handles.matfile.spikes(1),get(handles.tag_egm,'ylim'),'--r','linewidth',2)
end

if ~isequal(xx_lim,[0 1])
    set(handles.tag_egm,'xlim',[min(xx_lim) max(xx_lim)]);
    zoom off
end
if isnumeric(i_chan)&&~isempty(yy_lim)
    set(handles.tag_egm,'ylim',yy_lim);
end

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

handles.tag_map.BackgroundColor = [.8 .2 .2];

delete(findobj(handles.figure1,'tag','legend'))
delete(findobj(handles.tag_egm,'type','line'))
handles.cclines.handles = [];
handles.cclines.color = [];
cla(handles.tag_mesh)
cla(handles.tag_egm)
guidata(hObject,handles)



ib = str2double(get(handles.tag_beat_number,'string'));
str = get(handles.tag_menu,'string');
istr = get(handles.tag_menu,'value');
if isequal(str{istr},'Activ. T')
    if handles.tag_median_map.Value
        values = handles.matfile.Markers.dt - handles.matfile.spikes(:)*ones(1,size(handles.matfile.Markers.dt,2));
        values = nanmedian(values,1);
    else
        values = handles.matfile.Markers.dt(ib,:);
        values = values - handles.matfile.spikes(ib);
    end
    colorLabel = 'AT (ms)';
    set(handles.tag_slider,'value',0,'SliderStep',[0.01  0.1])
    handles.tag_modify_type.Value = 1;
end
if isequal(str{istr},'Repol. T')
    if handles.tag_median_map.Value
        values = handles.matfile.Markers.rt_Wyatt - handles.matfile.spikes(:)*ones(1,size(handles.matfile.Markers.dt,2));
        values = nanmedian(values,1);
    else
        values = handles.matfile.Markers.rt_Wyatt(ib,:);
        values = values -  handles.matfile.spikes(ib);
    end
    colorLabel = 'RT (ms)';
    set(handles.tag_slider,'value',0,'SliderStep',[0.01  0.1])
    handles.tag_modify_type.Value = 2;
end
if isequal(str{istr},'ARI')
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    if handles.tag_median_map.Value
        values = nanmedian(handles.matfile.Markers.rt_Wyatt -handles.matfile.Markers.dt,1);
    else
        values = handles.matfile.Markers.rt_Wyatt(ib,:) -handles.matfile.Markers.dt(ib,:);
    end
    colorLabel = 'ARI (ms)';
end


if isequal(str{istr},'Potential')
    handles.tag_median_map.Value = 0;
    t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);
    set(handles.tag_slider,'SliderStep',[1/range(t) 0.1])
    t0 = get(handles.tag_slider,'value')*t(end);
    [~,n] = min(abs(t-t0));
    values =  handles.matfile.signals(n,:)/max(abs(handles.matfile.signals(n,:)));
    colorLabel = ['Volt Norm (',num2str(t0,'%1.0f'),'ms)'];
end
if isequal(str{istr},'dV/dt')
    handles.tag_median_map.Value = 0;
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
        if handles.tag_median_map.Value
            values = handles.matfile.RVI.(name)(length(handles.matfile.spikes)+1,:);
        else
            values = handles.matfile.RVI.(name)(ib,:);
        end
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
    handles.tag_median_map.Value = 0;
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_AT(ib,:);
    colorLabel = 'Disp(AT) (ms)';
end
if isequal(str{istr},'DispRT')
    handles.tag_median_map.Value = 0;
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_RT(ib,:);
    colorLabel = 'Disp(RT) (ms)';
end
if isequal(str{istr},'DispARI')
    handles.tag_median_map.Value = 0;
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_ARI(ib,:);
    colorLabel = 'Disp(ARI) (ms)';
end
if isequal(str{istr},'DispAT(N)')
    handles.tag_median_map.Value = 0;
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_AT_norm(ib,:);
    colorLabel = 'Disp(AT) (n.u.)';
end
if isequal(str{istr},'DispRT(N)')
    handles.tag_median_map.Value = 0;
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_RT_norm(ib,:);
    colorLabel = 'Disp(RT) (n.u.)';
end
if isequal(str{istr},'DispARI(N)')
    handles.tag_median_map.Value = 0;
    set(handles.tag_slider,'value',0,'SliderStep',[0.01 0.1])
    values = handles.matfile.Disp_Local.range_ARI_norm(ib,:);
    colorLabel = 'Disp(ARI) (n.u.)';
end

if ~isnumeric(eventdata)||isempty(eventdata)
    Clim = [prctile(values,5) prctile(values,95)];
else
    Clim = eventdata; % from "reset"
end

if sum(isnan(Clim))>0
    warning('All nan')
    return
end


dcm_obj = datacursormode(handles.figure1);
set(dcm_obj,'UpdateFcn',@myfunctioncursor_markers_CARTO,'SnapToDataVertex','on')
cla(handles.figure1)
axes(handles.tag_mesh);

%% Meshes
MV.faces = handles.matfile.geo.Cmesh.triangles;
MV.vertices =  handles.matfile.geo.Cmesh.vertices;

if handles.tag_project.Value
    if ~isfield(handles.matfile.geo.Cmesh,'index_xyz');
        index_xyz = nan(size(handles.matfile.geo.Cmesh.vertices,1),1);
        for i = 1:size(handles.matfile.geo.Cmesh.vertices,1);
            [~,index_xyz(i)] = min(sum(abs( handles.matfile.geo.xyz - ones(size(handles.matfile.geo.xyz,1),1)*handles.matfile.geo.Cmesh.vertices(i,:) ),2));
        end
        handles.matfile.geo.Cmesh.index_xyz = index_xyz;
    end
    
    do_interpol = 1; % this may be brought to the GUI
    v = values(handles.matfile.geo.Cmesh.index_xyz(:));v = v(:);
    iiok = ~isnan(v);
    if do_interpol && sum(iiok)>0
        P = handles.matfile.geo.Cmesh.vertices;
        F = scatteredInterpolant(P(iiok,1),P(iiok,2),P(iiok,3),v(iiok));
        v_interp = F([P(:,1),P(:,2),P(:,3)]);
        MV.faceVertexCData = v_interp;
    else
        MV.faceVertexCData = v;
    end
else
    MV.faceVertexCData = 1;
end

% MV.faceVertexCData = values(:);
%
%% Modify the beginning of this function to just update the colors
hp = findobj(handles.tag_mesh,'type','patch','tag','Mesh_Ventricles');
if ~isempty(hp)
    hp.FaceVertexCData = MV.faceVertexCData;
else
    %% Ventricles
    meshObject = patch(MV,'parent',handles.tag_mesh);
    meshObject.Tag = 'Mesh_Ventricles';
    light('Position',[0 0 1],'Style','infinite')
    
    if exist([pwd,'\GUI_egm_mFiles\Geo_Chann\colormap_carto.mat'],'file')>0
        load([pwd,'\GUI_egm_mFiles\Geo_Chann\colormap_carto.mat'])
    else
        cmap = jet;
        cmap = cmap(end:-1:1,:);
    end
    
    if ~handles.tag_project.Value
        x = values;
        hold(handles.tag_mesh,'on');
        a2 = gobjects(size(x));
        for j = 1:length(x)
            a2(j) = plot3(handles.matfile.geo.xyz(j,1),handles.matfile.geo.xyz(j,2),handles.matfile.geo.xyz(j,3),'ok');
            xp = [x(j)-Clim(1)]/range(Clim);
            xp = round(xp*size(cmap,1)+1);
            xp(xp<1)=1;xp(xp>size(cmap,1))=size(cmap,1);
            if ~isnan(xp)
                a2(j).MarkerFaceColor = cmap(xp,:);
            end
        end
        axis off image
        set(a2,'MarkerSize',8,'MarkerEdgeColor','k');
        set(a2(isnan(x)),'visible','off');
        colormap(gca,cmap) ;
        set(gca,'clim',Clim)
        
        meshObject.FaceColor = [1 1 1]*.75;
        meshObject.EdgeColor = 'none';
        
        if handles.tag_CBreset_only.Value
            set(a2([x>Clim(2)|x<Clim(1)]),'visible','off');
        end
    else
        set(meshObject,'edgecolor','none','facecolor','interp','facelighting','gouraud',...
            'ambientstrength',0.9,'specularcolorreflectance',0.4,'specularexponent',3,'specularstrength',0.3);
        
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
    
    
    colormap(cmap);
    
    set(handles.tag_mesh,'clim',Clim)
    axis(handles.tag_mesh,'off','image')
    
    if get(handles.tag_nodes,'value')
        set(meshObject,'marker','o','MarkerFaceColor',[1 1 1]*.4,'markersize',4)
    else
        set(meshObject,'marker','none','MarkerFaceColor','none')
    end
    
    handles.tag_map.BackgroundColor = [1 1 1]*.94;
    
end


if handles.tag_menu.Value==6
    D = str2double(handles.tag_RVI_d.String);
    a = [xlim(handles.tag_mesh) ylim(handles.tag_mesh) zlim(handles.tag_mesh)];
    [x,y,z] = sphere;
    hold(handles.tag_mesh,'on');
    hs = surf(D*x+a(2)-[a(2)-a(1)]/6,D*y+a(4)-[a(4)-a(3)]/6,D*z+a(6)-[a(6)-a(5)]/6,'parent',handles.tag_mesh);
    clear x y z a
end

% -
if handles.tag_project.Value
    handles.tag_slider_alpha.Value = 1;
else
    handles.tag_slider_alpha.Value = .2;
end
tag_slider_alpha_Callback(hObject, [], handles)


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


if handles.tag_plot_filt_raw.Value
    S = handles.matfile.signals;
elseif handles.tag_plot_filt.Value
    S = handles.matfile.signals_proc;
elseif handles.tag_plot_filt_AT.Value
    S = handles.matfile.signals_proc_AT;
end

% Find beat
dd = pos(1) - handles.matfile.spikes;
dd(dd<0)=inf;
[~,ib] = min(dd);
ib0 = str2double(handles.tag_beat_number.String);

if handles.tag_modify_all_beat.Value
    samp_all = round([pos(1)+handles.matfile.spikes(1)-handles.matfile.spikes(ib)+handles.matfile.spikes]/1000*handles.matfile.ParamSig.frequency);
    samp0 = round([pos(1)+handles.matfile.spikes(1)]/1000*handles.matfile.ParamSig.frequency);
    ib_all = 1:length(handles.matfile.spikes);
else
    samp_all = round([pos(1)+handles.matfile.spikes(1)]/1000*handles.matfile.ParamSig.frequency);
    samp0 = samp_all;
    ib_all = ib
end
clear ib
% Find channel
ll = findobj(handles.figure1,'tag','legend');
ic_all = nan(size(ll.String));
for i=1:length(ll.String)
    ic_all(i) = str2double(ll.String{i}(6:end));
end
if handles.tag_modify_all_chan.Value
    ic = ic_all;
else
    ic = ic_all(S(samp0,ic_all)==pos(2));
end
clear samp0
% =

if isempty(ic)
    warning('Channel not identified')
    return
end
T = round(str2double(get(handles.tag_window,'string'))/1000*handles.matfile.ParamSig.frequency/2);
t = [1:size(handles.matfile.signals,1)]/handles.matfile.ParamSig.frequency*1000 - handles.matfile.spikes(1);

for i_samp = 1:length(samp_all);
    samp = samp_all(i_samp);
    ib = ib_all(i_samp);
    
    Ws = samp + [-T : T]; % modified 29/11/2016
    Ws(Ws<1)=[];
    X = S(Ws,ic);
    if get(handles.tag_modify_type,'value')==1
        
        % Modify activation
        Xd = diff(X);Xd2 = diff(Xd);Xd3 = diff(Xd2);
        ii = Xd2(1:end-1,:).*Xd2([2 : end],:)<0 & Xd([1: end-2],:)<0 & Xd3([1: end],:)>0; % zeros of II derivative & I derivative negative & III derivative > 0
        
        for j = 1:size(ii,2)
            tt = find(ii(:,j));
            if ~isempty(tt)
                [~,kk]=min(Xd(tt,j));
                handles.matfile.Markers.dt(ib,ic(j)) = [tt(kk)+Ws(1)]/handles.matfile.ParamSig.frequency*1000;
                % plot new marker
                dt = handles.matfile.Markers.dt(ib,ic(j)) - handles.matfile.spikes(1);
                % AT
                dtp = round(handles.matfile.Markers.dt(ib,ic(j))/1000*handles.matfile.ParamSig.frequency);
                plot(handles.tag_egm,t(dtp),S(dtp,ic(j)),'xr','markersize',12)
            end
        end
        %         plot(handles.tag_egm,pos(1)-str2double(get(handles.tag_window,'string'))/2*[1 1],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
        %         plot(handles.tag_egm,pos(1)+str2double(get(handles.tag_window,'string'))/2*[1 1],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
        plot(handles.tag_egm,[t(Ws(1)) t(Ws(1))],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
        plot(handles.tag_egm,[t(Ws(end)) t(Ws(end))],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
        
        clear x xd*
    else
        
        Xd = diff(X);Xd2 = diff(Xd);Xd3 = diff(Xd2);
        ii = Xd2([1:end-1],:).*Xd2([2 : end],:)<0 & Xd([1: end-2],:)>0 & Xd3(1:end,:)<=0;
        
        for j = 1:size(ii,2)
            tt = find(ii(:,j));
            
            if ~isempty(tt)
                [~,kk]=max(Xd(tt,j));
                handles.matfile.Markers.rt_Wyatt(ib,ic(j)) = [tt(kk)+Ws(1)]/handles.matfile.ParamSig.frequency*1000;
                
                % plot new marker
                rt = handles.matfile.Markers.rt_Wyatt(ib,ic(j)) - handles.matfile.spikes(1);
                % RT
                rtp = round(handles.matfile.Markers.rt_Wyatt(ib,ic(j))/1000*handles.matfile.ParamSig.frequency)+1;
                plot(handles.tag_egm,t(rtp),S(rtp,ic(j)),'xr','markersize',12)
            end
        end
        %         plot(handles.tag_egm,pos(1)-str2double(get(handles.tag_window,'string'))/2*[1 1],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
        %         plot(handles.tag_egm,pos(1)+str2double(get(handles.tag_window,'string'))/2*[1 1],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
        plot(handles.tag_egm,[t(Ws(1)) t(Ws(1))],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
        plot(handles.tag_egm,[t(Ws(end)) t(Ws(end))],range(get(handles.tag_egm,'ylim'))/5*[-1 1]+pos(2),'--k')
        
        clear x xd*
    end
end


if handles.tag_project.Value
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
else
    a2 = findobj(handles.tag_mesh,'marker','o','markersize',8);
    a2 = a2(end:-1:1);
    
    Clim = get(handles.tag_mesh,'clim');
    cmap = colormap;
    if handles.tag_modify_type.Value==1
        x = handles.matfile.Markers.dt - handles.matfile.spikes(:)*ones(1,size(handles.matfile.Markers.dt,2));
    else
        x = handles.matfile.Markers.rt_Wyatt - handles.matfile.spikes(:)*ones(1,size(handles.matfile.Markers.dt,2));
    end
    
    for j = 1:length(ic)
        xp = [x(ib0,ic(j))-Clim(1)]/range(Clim);
        xp = round(xp*size(cmap,1)+1);
        xp(xp<1)=1;xp(xp>size(cmap,1))=size(cmap,1);
        if ~isnan(xp)
            a2(ic(j)).MarkerFaceColor = cmap(xp,:);
        else
            %                 a2(ic(j)).MarkerFaceColor = 'k';
            %                 a2(ic(j)).Marker = 'x';
            %             a2(ic(j)).MarkerSize = 3;
            set(a2(ic(j)),'visible','off');
        end
    end
    
    
end


guidata(hObject,handles)
delete(findall(gcf,'Type','hggroup'));



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

po = findobj(handles.tag_mesh,'type','patch');
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
if handles.tag_median_map.Value
    ib = 1:length(handles.matfile.spikes);
else
    ib = str2double(handles.tag_beat_number.String);
end

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
    
    clear ib
    
elseif isequal(get(info_struct.Target.Parent,'type'),'axes')
    if handles.tag_plot_filt_raw.Value
        S = handles.matfile.signals;
    elseif handles.tag_plot_filt.Value
        S = handles.matfile.signals_proc;
    elseif handles.tag_plot_filt_AT.Value
        S = handles.matfile.signals_proc_AT;
    end
    
    pos = info_struct.Position;
    dt =pos(1)-handles.matfile.spikes;
    dt(dt<0)=inf;
    [~,ib] = min(dt);
    clear dt
    
    % RT
    x = handles.matfile.Markers.rt_Wyatt(ib,:)-handles.matfile.spikes(1);
    ii = find(x==pos(1));
    jj = zeros(1,length(ii));
    for i = 1:length(ii)
        jj(i) = (S(round(handles.matfile.Markers.rt_Wyatt(ib,ii(i))/1000*handles.matfile.ParamSig.frequency),ii(i))==pos(2));
    end
    ic = ii(logical(jj));
    
    if ~isempty(ic)
        if handles.tag_project.Value==0
            po = findobj(handles.tag_mesh,'marker','o','markersize',8);
            po = po(end:-1:1);
            po(ic).Marker = 'x';
        else
            po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
            po.FaceVertexCData = handles.matfile.Markers.rt_Wyatt(ib,:)'-handles.matfile.spikes(ib);
        end
        handles.tag_menu.Value = 2;
        
        if handles.tag_modify_all_beat.Value
            ib = 1:length(handles.matfile.spikes);
        end
        if handles.tag_modify_all_chan.Value
            pp = findobj(handles.figure1,'tag','legend');
            ic = nan(size(pp.String));
            for i = 1:length(pp.String)
                ic(i) = str2double(pp.String{i}(6:end));
            end
            clear i
        end
        for i = 1:length(ic)
            handles.matfile.Markers.rt_Wyatt(ib,ic(i)) = nan;
        end
    end
    clear ic;
    
    % AT
    x = handles.matfile.Markers.dt(ib,:)-handles.matfile.spikes(1);
    ii = find(x==pos(1));
    jj = zeros(1,length(ii));
    for i = 1:length(ii)
        jj(i) = (S(round(handles.matfile.Markers.dt(ib,ii(i))/1000*handles.matfile.ParamSig.frequency),ii(i))==pos(2));
    end
    ic = ii(logical(jj));
    
    
    if ~isempty(ic)
        if handles.tag_project.Value==0
            po = findobj(handles.tag_mesh,'marker','o','markersize',8);
            po = po(end:-1:1);
            po(ic).Marker = 'x';
        else
            po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
            po.FaceVertexCData = handles.matfile.Markers.dt(ib,:)'-handles.matfile.spikes(ib);
        end
        handles.tag_menu.Value = 1;
        %
        if handles.tag_modify_all_beat.Value
            ib = 1:length(handles.matfile.spikes);
        end
        if handles.tag_modify_all_chan.Value
            pp = findobj(handles.figure1,'tag','legend');
            ic = nan(size(pp.String));
            for i = 1:length(pp.String)
                ic(i) = str2double(pp.String{i}(6:end));
            end
            clear i
        end
        for i = 1:length(ic)
            handles.matfile.Markers.dt(ib,ic(i)) = nan;
        end
    end
    clear ic;
end
guidata(hObject,handles)

iv = handles.tag_hold_on.Value;
handles.tag_hold_on.Value=0;

tag_plot_egm_Callback(hObject,handles.tag_find_chan.Value,handles);
handles.tag_hold_on.Value = iv;


% --- Executes on button press in tag_modify_all_chan.
function tag_modify_all_chan_Callback(hObject, eventdata, handles)
% hObject    handle to tag_modify_all_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_modify_all_chan



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

pa = findobj(handles.tag_egm,'type','line');
if ~isempty(pa)
    if ib == length(handles.matfile.spikes)
        xlim(handles.tag_egm,[handles.matfile.spikes(ib)-500 handles.matfile.spikes(ib)+1500]-handles.matfile.spikes(1));
    else
        xlim(handles.tag_egm,[handles.matfile.spikes(ib)-500 handles.matfile.spikes(ib+1)+500]-handles.matfile.spikes(1));
    end
end
plot(handles.tag_egm,[handles.matfile.spikes(ib)]*[1 1]-handles.matfile.spikes(1),get(handles.tag_egm,'ylim'),'--r','linewidth',2)



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




% --- Executes on button press in tag_RVI.
function tag_RVI_Callback(hObject, eventdata, handles)
% hObject    handle to tag_RVI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% RVI and Dispersion calculated using the median markers are reported in
% length(spikes)+1 row.

if ~isfield(handles.matfile,'Disp_Local')
    handles.matfile.Disp_Local.range_AT = nan(length(handles.matfile.spikes)+1,size(handles.matfile.Markers.dt,2));
    handles.matfile.Disp_Local.range_ARI = nan(length(handles.matfile.spikes)+1,size(handles.matfile.Markers.dt,2));
    handles.matfile.Disp_Local.range_RT =nan(length(handles.matfile.spikes)+1,size(handles.matfile.Markers.dt,2));
    handles.matfile.Disp_Local.range_AT_norm = nan(length(handles.matfile.spikes)+1,size(handles.matfile.Markers.dt,2));
    handles.matfile.Disp_Local.range_ARI_norm = nan(length(handles.matfile.spikes)+1,size(handles.matfile.Markers.dt,2));
    handles.matfile.Disp_Local.range_RT_norm = nan(length(handles.matfile.spikes)+1,size(handles.matfile.Markers.dt,2));
end


if handles.tag_median_map.Value
    dt = nanmedian(handles.matfile.Markers.dt - handles.matfile.spikes(:)*ones(1,size(handles.matfile.Markers.dt,2)),1);
    rt = nanmedian(handles.matfile.Markers.rt_Wyatt - handles.matfile.spikes(:)*ones(1,size(handles.matfile.Markers.dt,2)),1);
    ib = length(handles.matfile.spikes)+1; % write RVI and DISP values from median markers here
else
    ib = str2double(handles.tag_beat_number.String);
    dt = handles.matfile.Markers.dt(ib,:) - handles.matfile.spikes(ib);
    rt = handles.matfile.Markers.rt_Wyatt(ib,:) - handles.matfile.spikes(ib);
end

D = str2double(handles.tag_RVI_d.String);
rvi_type = handles.tag_RVI_type.Value;
xyz = handles.matfile.geo.xyz;
values = nan(length(dt),1);
Nnodes = nan(size(values));


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
        handles.matfile.RVI.(['min_D_',handles.tag_RVI_d.String]) = nan(length(handles.matfile.spikes)+1,size(handles.matfile.signals,2));
        handles.matfile.RVI.(['min_D_',handles.tag_RVI_d.String,'_Nnodes']) = nan(length(handles.matfile.spikes)+1,size(handles.matfile.signals,2));
    end
    handles.matfile.RVI.(['min_D_',handles.tag_RVI_d.String])(ib,:) = values(:).';
    handles.matfile.RVI.(['min_D_',handles.tag_RVI_d.String,'_Nnodes'])(ib,:) = Nnodes(:).';
elseif rvi_type==2
    if ~isfield(handles.matfile.RVI,['avg_D_',handles.tag_RVI_d.String]);
        handles.matfile.RVI.(['avg_D_',handles.tag_RVI_d.String]) = nan(length(handles.matfile.spikes)+1,size(handles.matfile.signals,2));
        handles.matfile.RVI.(['avg_D_',handles.tag_RVI_d.String,'_Nnodes']) = nan(length(handles.matfile.spikes)+1,size(handles.matfile.signals,2));
    end
    handles.matfile.RVI.(['avg_D_',handles.tag_RVI_d.String])(ib,:) =values(:).';
    handles.matfile.RVI.(['avg_D_',handles.tag_RVI_d.String,'_Nnodes'])(ib,:) = Nnodes(:).';
end


handles.tag_menu.Value = 6;
guidata(hObject,handles);
%
tag_map_Callback(hObject,[],handles);

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
aa = plot3(handles.matfile.geo.xyz(ic,1),handles.matfile.geo.xyz(ic,2),handles.matfile.geo.xyz(ic,3),'x');
set(aa,'markerfacecolor',[1 1 1]*.6,'markeredgecolor','w','markersize',20,'linewidth',2)

%
tag_plot_egm_Callback(hObject,ic,handles);



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
% pcbar = findobj(handles.tag_mesh,'tag','colorbar');

if handles.tag_project.Value
    set(handles.tag_mesh,'Clim',Clim);
    cbar = findobj(handles.figure1,'tag','Colorbar');
    set(cbar,'xlim',Clim)
else
    tag_map_Callback(hObject,Clim, handles);
end

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


if handles.tag_export_fig.Value==2 % RVI+VT
    
    tag_fig_summary_Callback(hObject, eventdata, handles)
    handles.tag_export_fig.Value=1;
    
elseif handles.tag_export_fig.Value==3 % Just export the current map
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
    
elseif handles.tag_export_fig.Value==4 % Just export the current map
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


% --- Executes on button press in tag_segment.
function tag_segment_Callback(hObject, eventdata, handles)
% hObject    handle to tag_segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_segment_mesh(handles);


% --- Executes on button press in tag_project.
function tag_project_Callback(hObject, eventdata, handles)
% hObject    handle to tag_project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_project
tag_map_Callback(hObject,[], handles);


% --- Executes on button press in tag_median_map.
function tag_median_map_Callback(hObject, eventdata, handles)
% hObject    handle to tag_median_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_median_map


% --- Executes on button press in tag_plot_filt_AT.
function tag_plot_filt_AT_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_filt_AT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_filt_AT
if handles.tag_plot_filt_AT.Value
    handles.tag_plot_filt_raw.Value = 0;
    handles.tag_plot_filt.Value = 0;
    tag_plot_egm_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in tag_plot_filt_raw.
function tag_plot_filt_raw_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_filt_raw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_filt_raw
% --- Executes on button press in tag_plot_filt.
if handles.tag_plot_filt_raw.Value
    handles.tag_plot_filt_AT.Value = 0;
    handles.tag_plot_filt.Value = 0;
    tag_plot_egm_Callback(hObject, eventdata, handles)
end

function tag_plot_filt_Callback(hObject, eventdata, handles)
% hObject    handle to tag_plot_filt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_plot_filt
if handles.tag_plot_filt.Value
    handles.tag_plot_filt_AT.Value = 0;
    handles.tag_plot_filt_raw.Value = 0;
    tag_plot_egm_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in tag_modify_all_beat.
function tag_modify_all_beat_Callback(hObject, eventdata, handles)
% hObject    handle to tag_modify_all_beat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_modify_all_beat


% --- Executes on slider movement.
function tag_slider_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to tag_slider_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

pp = findobj(handles.tag_mesh,'tag','Mesh_Ventricles');
pp.FaceAlpha = handles.tag_slider_alpha.Value;

% --- Executes during object creation, after setting all properties.
function tag_slider_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_slider_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in tag_CBreset_only.
function tag_CBreset_only_Callback(hObject, eventdata, handles)
% hObject    handle to tag_CBreset_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_CBreset_only


% --------------------------------------------------------------------
function tag_menu_tool_Callback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function tag_import_point_Callback(hObject, eventdata, handles)
% hObject    handle to tag_import_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FILENAME, PATHNAME] = uigetfile('*txt','MultiSelect', 'on');
Pos = cell(1,length(FILENAME));
if ~iscell(FILENAME)
    FILENAME = {FILENAME};
end
for i=1:length(FILENAME)
    loadname_pos = [PATHNAME,FILENAME{i}];
    
    
    
    if ~isempty(strfind(loadname_pos,'NAVISTAR_CONNECTOR_Eleclectrode_Positions.txt'));
        ii1 = findstr(FILENAME{i},'_P');ii1 = ii1(1);
        ii2 = findstr(FILENAME{i},'_');ii2(ii2<=ii1)=[];ii2 = ii2(1);
            pnumber = FILENAME{i}(ii1+1:ii2-1);
        
        ll = {'M1','M2','M3','M4'};
        chan =  listdlg('liststring',ll,'PromptString',['Select ELECTRODES for P=',pnumber],'SelectionMode','multiple');
        Pos{i} = LoadElPositionFromCarto(loadname_pos);
        xyz = squeeze(nanmean(Pos{i}.xyz(:,:,chan)))';
        
    elseif ~isempty(strfind(loadname_pos,'_20_POLE_A_CONNECTOR_Eleclectrode_Positions.txt'))|~isempty(strfind(loadname_pos,'_20_POLE_B_CONNECTOR_Eleclectrode_Positions.txt'))
        ii1 = findstr(FILENAME{i},'_P');ii1 = ii1(1);
        ii2 = findstr(FILENAME{i},'_');ii2(ii2<=ii1)=[];ii2 = ii2(1);
            pnumber = FILENAME{i}(ii1+1:ii2-1);
        ll = {'Penta-1','Penta-2','Penta-3','Penta-4','Penta-5','Penta-6','Penta-7','Penta-8','Penta-9','Penta-10',...
            'Penta-11','Penta-12','Penta-13','Penta-14','Penta-15','Penta-16','Penta-17','Penta-18','Penta-19','Penta-20'};
        chan =  listdlg('liststring',ll,'PromptString',['Select ELECTRODES for P=',pnumber],'SelectionMode','multiple');
        Pos{i} = LoadElPositionFromCarto(loadname_pos);
        xyz = squeeze(nanmean(Pos{i}.xyz(:,:,chan)))';
    elseif ~isempty(strfind(loadname_pos,'Sites.txt'))
        fileID = fopen(loadname_pos,'r');
        l = fgetl(fileID);
        dataArray = textscan(fileID,repmat('%f',[1 14]),'CollectOutput',1);
        fopen(fileID);
        ll = [num2str(dataArray{1}(:,[1])),repmat('-',[size(dataArray{1},1),1]),num2str(dataArray{1}(:,[3]))];
        chan =  listdlg('liststring',ll,'PromptString',['Select ABLATION point'],'SelectionMode','multiple');
        xyz = dataArray{1}(chan,4:6);
    elseif ~isempty(strfind(loadname_pos,'Positions_OnAnnotation.txt'));
        ii1 = findstr(FILENAME{i},'_P');ii1 = ii1(1);
        ii2 = findstr(FILENAME{i},'_');ii2(ii2<=ii1)=[];ii2 = ii2(1);
        pnumber = FILENAME{i}(ii1+1:ii2-1);
        
        Pos{i} = LoadElPositionFromCarto(loadname_pos);
        N= size(Pos{i}.xyz,3);
        ll = cell(1,N);
        for k = 1:N
        ll{k} = ['Chan-',num2str(k)];
        end
           
        chan =  listdlg('liststring',ll,'PromptString',['Select ELECTRODES for P=',pnumber],'SelectionMode','multiple');
        
        if length(chan)>1
        xyz = squeeze(Pos{i}.xyz(:,:,chan))';
        else
        xyz = squeeze(Pos{i}.xyz(:,:,chan));    
        end
    else
        error('Point should be either a mapping or a penta-array position point');
    end
    
    
    hold(handles.tag_mesh,'on');
%     plot3(handles.tag_mesh,Pos{i}.xyz_OnAnnot(chan,1),Pos{i}.xyz_OnAnnot(chan,2),Pos{i}.xyz_OnAnnot(chan,3),'x','markersize',20,'linewidth',2,'color','k');
    plot3(handles.tag_mesh,xyz(:,1),xyz(:,2),xyz(:,3),'x','markersize',20,'linewidth',2,'color','k');

end


% --------------------------------------------------------------------
function tag_measure_dist_Callback(hObject, eventdata, handles)
% hObject    handle to tag_measure_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_measure_distance(handles.figure1)


% --------------------------------------------------------------------
function tag_export_map_fig_Callback(hObject, eventdata, handles)
% hObject    handle to tag_export_map_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


hfn = figure;

cb = findobj(handles.figure1,'type','colorbar');

ax = copyobj([handles.tag_mesh,cb],hfn);


cc = findobj(ax,'type','surface');

delete(cc)
axis(ax(1),'image');
set(ax(1),'position',[.1 .2 .85 .85]);
set(ax(2),'position',[.21 .05 .6 .05],'fontsize',11);

[n,p] = uiputfile('.fig','Save As',handles.matfile.ParamSig.name_file(1:end-4));
saveas(hfn,[p,n],'fig');
