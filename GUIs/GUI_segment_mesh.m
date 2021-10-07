function varargout = GUI_segment_mesh(varargin)
% GUI_SEGMENT_MESH MATLAB code for GUI_segment_mesh.fig
%      GUI_SEGMENT_MESH, by itself, creates a new GUI_SEGMENT_MESH or raises the existing
%      singleton*.
%
%      H = GUI_SEGMENT_MESH returns the handle to a new GUI_SEGMENT_MESH or the handle to
%      the existing singleton*.
%
%      GUI_SEGMENT_MESH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SEGMENT_MESH.M with the given input arguments.
%
%      GUI_SEGMENT_MESH('Property','Value',...) creates a new GUI_SEGMENT_MESH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_segment_mesh_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_segment_mesh_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_segment_mesh

% Last Modified by GUIDE v2.5 03-Oct-2017 11:08:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_segment_mesh_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_segment_mesh_OutputFcn, ...
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


% --- Executes just before GUI_segment_mesh is made visible.
function GUI_segment_mesh_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_segment_mesh (see VARARGIN)

% Choose default command line output for GUI_segment_mesh
handles.output = hObject;

handles.handles_ori = varargin{1};

tag_mesh = copyobj(handles.handles_ori.tag_mesh,handles.figure1,'legacy');
handles.tag_mesh = tag_mesh;
cc =colormap(handles.handles_ori.tag_mesh);
colormap(tag_mesh,cc);
set(handles.tag_mesh,'position',[.05 .05 .7 .8])
cc = findobj(handles.tag_mesh,'type','patch','tag','Mesh_Ventricles');
set(cc,'edgecolor','k','marker','o','markeredgecolor','k','markerfacecolor','flat','facecolor',[1 1 1]*.7,'markerfacecolor','w','markersize',6);

% list
label_list = {'valve_plane','basal_anterior_lv','mid_anterior_lv','basal_anterolateral_lv',...,
    'mid_anterolateral_lv','basal_inferolateral_lv','mid_inferolateral_lv','basal_inferior_lv',...,
    'mid_inferior_lv','apical_anterior_lv','apical_inferior_lv', 'basal_anterior_rv','mid_anterior_rv',...,
    'basal_anterolateral_rv','mid_anterolateral_rv','basal_inferolateral_rv', 'mid_inferolateral_rv',...,
    'basal_inferior_rv','mid_inferior_rv','apical_anterior_rv', 'apical_inferior_rv','LV_ALL','RV_ALL'};
set(handles.tag_list_label,'string',label_list);


handles.tag_n_nodes.String = num2str(size(handles.handles_ori.matfile.geo.Nodes_ventricles,1));

if isfield(handles.handles_ori.matfile.geo,'Labels_anatomy')
    handles.Labels_anatomy = handles.handles_ori.matfile.geo.Labels_anatomy;
    vv = fieldnames(handles.Labels_anatomy);
    [~,iiA,iiB] =  intersect(handles.tag_list_label.String,vv,'stable');
  
    for i = 1:length(iiA)
        handles.tag_list_label.String{iiA(i)} = ['<HTML><FONT color="red">',handles.tag_list_label.String{iiA(i)},' [',num2str(length(handles.Labels_anatomy.(vv{iiB(i)}).index)),']</Font></html>'];
    end
    handles.tag_list_label.Value = iiA;
    handles.tag_show_roi.Value =1;
    tag_show_roi_Callback(hObject,[], handles);
end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_segment_mesh wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_segment_mesh_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in tag_add_point.
function tag_add_point_Callback(hObject, eventdata, handles)
% hObject    handle to tag_add_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_add_point
if handles.tag_add_point.Value;
    datacursormode on;
    dcm_obj = datacursormode(handles.figure1);
    set(dcm_obj,'UpdateFcn',@myfunctioncursor_segment,'SnapToDataVertex','on');
else
    datacursormode off
end




% --- Executes on button press in tag_grow.
% function tag_grow_Callback(hObject, eventdata, handles)
% % hObject    handle to tag_grow (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % keyboard
% cc = findobj(handles.tag_mesh,'tag','segment');
%
% xyz0 = [[cc.XData]' [cc.YData]' [cc.ZData]'];
% [~,ii] = unique(sum(xyz0,2));
% xyz = xyz0(ii,:);
%
% iiko = setdiff([1:size(xyz0,1)],ii);
% delete(cc(iiko));
%
% xyz0 = handles.handles_ori.matfile.geo.Nodes_ventricles;
% iipoints = xyz0(:,1)<=max(xyz(:,1)) & xyz0(:,1)>=min(xyz(:,1)) & xyz0(:,2)<=max(xyz(:,2)) & xyz0(:,2)>=min(xyz(:,2)) & xyz0(:,3)<=max(xyz(:,3)) & xyz0(:,3)>=min(xyz(:,3));
%
% iipoints = find(iipoints);
% for i = 1:length(iipoints)
% aa(i) = plot3(handles.tag_mesh,xyz0(iipoints(i),1),xyz0(iipoints(i),2),xyz0(iipoints(i),3),'or','markerfacecolor','g','markersize',8);
% aa(i).Tag = 'segment';
% end

% --- Executes on button press in tag_nodes.
function tag_nodes_Callback(hObject, eventdata, handles)
% hObject    handle to tag_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_nodes
po = findobj(handles.tag_mesh,'type','patch','FaceColor','interp');
if get(handles.tag_nodes,'value')
    set(po,'marker','o','MarkerFaceColor',[1 1 1]*.4,'markersize',4)
else
    set(po,'marker','none','MarkerFaceColor','none')
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

% --- Executes on button press in tag_mesh_atria.
function tag_mesh_atria_Callback(hObject, eventdata, handles)
% hObject    handle to tag_mesh_atria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_mesh_atria

if ~handles.tag_mesh_atria.Value;
    iic = findobj(handles.tag_mesh,'type','patch');
    iin_tag = {iic.Tag};
    [~,iiok] = setdiff(iin_tag,'Mesh_Ventricles','stable');
    set(iic(iiok),'visible','off')
else
    iic = findobj(handles.tag_mesh,'type','patch');
    iin_tag = {iic.Tag};
    [~,iiok] = setdiff(iin_tag,'Mesh_Ventricles','stable');
    set(iic(iiok),'visible','on')
end


% --- Executes on selection change in tag_list_label.
function tag_list_label_Callback(hObject, eventdata, handles)
% hObject    handle to tag_list_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_list_label contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_list_label


% --- Executes during object creation, after setting all properties.
function tag_list_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_list_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_assign.
function tag_assign_Callback(hObject, eventdata, handles)
% hObject    handle to tag_assign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cc = findobj(handles.tag_mesh,'tag','segment');
if isempty(cc)
   return 
end
lab_name = handles.tag_list_label.String{handles.tag_list_label.Value};
if isequal(lab_name(1:6),'<HTML>');
    ii1 = find(lab_name=='>');
    ii2 = find(lab_name=='[');
    lab_name = lab_name(ii1(2)+1:ii2(1)-2);
end
xyz0 = [[cc.XData]' [cc.YData]' [cc.ZData]'];
[~,ii] = intersect(sum(handles.handles_ori.matfile.geo.Nodes_ventricles,2),sum(xyz0,2));

% check vertices already assigned
if handles.tag_overwrite.Value==0
    if isfield(handles,'Labels_anatomy')
        vv = fieldnames(handles.Labels_anatomy);
        ii_existing = [];
        for j = 1:length(vv)
            ii_existing = [ii_existing;handles.Labels_anatomy.(vv{j}).index];
        end
        [ii2] = setdiff(ii,ii_existing);
    else
        ii2 = ii;
    end
    
    handles.Labels_anatomy.(lab_name).index = ii2;
    handles.Labels_anatomy.(lab_name).xyz = handles.handles_ori.matfile.geo.Nodes_ventricles(ii2,:);
else
    if isfield(handles,'Labels_anatomy')
        vv = fieldnames(handles.Labels_anatomy);
        if ~contains(lab_name,'ALL')
        for j = 1:length(vv)
            [~,iiA,iiB] =  intersect(handles.Labels_anatomy.(vv{j}).index,ii,'stable');
            if ~isempty(iiA)
                handles.Labels_anatomy.(vv{j}).index(iiA) = [];
                handles.Labels_anatomy.(vv{j}).xyz(iiA,:) = [];
            end
        end
        end
        handles.Labels_anatomy.(lab_name).index = ii;
        handles.Labels_anatomy.(lab_name).xyz = handles.handles_ori.matfile.geo.Nodes_ventricles(ii,:);
    end 
end

iistr = handles.tag_list_label.Value;
str = handles.tag_list_label.String;
str{iistr} = ['<HTML><FONT color="red">',lab_name,' [',num2str(length(ii)),']</Font></html>'];
handles.tag_list_label.String = str;

delete(cc);
plot3(handles.Labels_anatomy.(lab_name).xyz(:,1),handles.Labels_anatomy.(lab_name).xyz(:,2),handles.Labels_anatomy.(lab_name).xyz(:,3),'diamondk','markerfacecolor','k');

handles.tag_list_label.Value = 1:length(handles.tag_list_label.String);
tag_show_roi_Callback(hObject, [], handles);

% updating number of nodes
ii_existing = [];
vv = fieldnames(handles.Labels_anatomy);
for j = 1:length(vv)
    ii_existing = [ii_existing;handles.Labels_anatomy.(vv{j}).index];
end
nn = setdiff([1:size(handles.handles_ori.matfile.geo.Nodes_ventricles,1)],ii_existing);      
handles.tag_n_nodes.String = num2str(length(nn));

%%
guidata(hObject, handles);


% --- Executes on button press in tag_remove_all.
function tag_remove_all_Callback(hObject, eventdata, handles)
% hObject    handle to tag_remove_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cc = findobj(handles.tag_mesh,'tag','segment');
delete(cc);
cc = findobj(handles.tag_mesh,'markerfacecolor','k');
delete(cc);
cc = findobj(handles.tag_mesh,'markerfacecolor','g');
delete(cc);



% --- Executes on selection change in tag_menu.
function tag_menu_Callback(hObject, eventdata, handles)
% hObject    handle to tag_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tag_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tag_menu
if handles.tag_menu.Value==1 % Add
    cc = findobj(handles.tag_mesh,'tag','segment');
    
    xyz0 = [[cc.XData]' [cc.YData]' [cc.ZData]'];
    [~,ii] = unique(sum(xyz0,2));
    iiko = setdiff([1:size(xyz0,1)],ii);
    delete(cc(iiko));
    
    for i = 1:size(xyz0,1)
        aa(i) = plot3(handles.tag_mesh,xyz0(i,1),xyz0(i,2),xyz0(i,3),'or','markerfacecolor','g','markersize',8);
        aa(i).Tag = 'segment';
    end
    
    cc = findobj(handles.tag_mesh,'tag','segment');
    xyz0 = [[cc.XData]' [cc.YData]' [cc.ZData]'];
    [~,ii] = unique(sum(xyz0,2));
    iiko = setdiff([1:size(xyz0,1)],ii);
    delete(cc(iiko));
    
    ccdel = findobj(handles.tag_mesh,'markerfacecolor','k');
    delete(ccdel)
    delete(findall(handles.figure1,'Type','hggroup'));
    
end


if handles.tag_menu.Value==2 % Grow
    cc = findobj(handles.tag_mesh,'tag','segment');
    
    xyz0 = [[cc.XData]' [cc.YData]' [cc.ZData]'];
    [~,ii] = unique(sum(xyz0,2));
    xyz = xyz0(ii,:);
    
    iiko = setdiff([1:size(xyz0,1)],ii);
    delete(cc(iiko));
    
    xyz0 = handles.handles_ori.matfile.geo.Nodes_ventricles;
    iipoints = xyz0(:,1)<=max(xyz(:,1)) & xyz0(:,1)>=min(xyz(:,1)) & xyz0(:,2)<=max(xyz(:,2)) & xyz0(:,2)>=min(xyz(:,2)) & xyz0(:,3)<=max(xyz(:,3)) & xyz0(:,3)>=min(xyz(:,3));
    
    iipoints = find(iipoints);
    for i = 1:length(iipoints)
        aa(i) = plot3(handles.tag_mesh,xyz0(iipoints(i),1),xyz0(iipoints(i),2),xyz0(iipoints(i),3),'or','markerfacecolor','g','markersize',8);
        aa(i).Tag = 'segment';
    end
    
    cc = findobj(handles.tag_mesh,'tag','segment');
    
    
    xyz0 = [[cc.XData]' [cc.YData]' [cc.ZData]'];
    [~,ii] = unique(sum(xyz0,2));
    iiko = setdiff([1:size(xyz0,1)],ii);
    delete(cc(iiko));
    delete(findall(handles.figure1,'Type','hggroup'));
    
end

if handles.tag_menu.Value==3 % Remove
    datacursormode off;
    handles.tag_add_point.Value=0;
    
    cc = findobj(handles.tag_mesh,'markerfacecolor','g');
    % delete copies
    xyz0 = [[cc.XData]' [cc.YData]' [cc.ZData]'];
    [~,ii] = unique(sum(xyz0,2));
    iiko = setdiff([1:size(xyz0,1)],ii);
    delete(cc(iiko));
    
    ccdel = findobj(handles.tag_mesh,'markerfacecolor','k');
    
    xyz0 = [[cc.XData]' [cc.YData]' [cc.ZData]'];
    if ~isempty(ccdel)
        xyzDel = [[ccdel.XData]' [ccdel.YData]' [ccdel.ZData]'];
    else
        delete(findall(handles.figure1,'Type','hggroup'));
        return
    end
    
    [~,ii] = intersect(sum(xyz0,2),sum(xyzDel,2),'stable');
    delete(ccdel);
    
    delete(cc(ii));
    delete(findall(handles.figure1,'Type','hggroup'));
end

if handles.tag_menu.Value==4 % Add [min]
    cc = findobj(handles.tag_mesh,'tag','segment');
    
    xyz0 = [[cc.XData]' [cc.YData]' [cc.ZData]'];
    [~,ii] = unique(sum(xyz0,2));
    xyz = xyz0(ii,:);
    
    iiko = setdiff([1:size(xyz0,1)],ii);
    delete(cc(iiko));
    
    xyz0 = handles.handles_ori.matfile.geo.Nodes_ventricles;
    
    [~,im] = min(sqrt(sum((xyz0 - mean(xyz)).^2,2)));
    xyzC = xyz0(im,:);
    R = mean(sqrt(sum([xyz - xyzC].^2,2)));
    d = sqrt(sum([xyz0 - xyzC].^2,2));
    iipoints = d<R;
        
    iipoints = find(iipoints);
    for i = 1:length(iipoints)
        aa(i) = plot3(handles.tag_mesh,xyz0(iipoints(i),1),xyz0(iipoints(i),2),xyz0(iipoints(i),3),'or','markerfacecolor','g','markersize',8);
        aa(i).Tag = 'segment';
    end
    
    cc = findobj(handles.tag_mesh,'tag','segment');
    
    
    xyz0 = [[cc.XData]' [cc.YData]' [cc.ZData]'];
    [~,ii] = unique(sum(xyz0,2));
    iiko = setdiff([1:size(xyz0,1)],ii);
    delete(cc(iiko));
    delete(findall(handles.figure1,'Type','hggroup'));
    
end


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


% --- Executes on button press in tag_show_roi.
function tag_show_roi_Callback(hObject, eventdata, handles)
% hObject    handle to tag_show_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.tag_show_roi.Value;
    cmap = lines;
    map = zeros(size(handles.handles_ori.matfile.geo.Nodes_ventricles,1),1);
    lab_name_tot = handles.tag_list_label.String(handles.tag_list_label.Value);
    iiko = cellfun(@isempty,strfind(lab_name_tot,'HTML'));
    lab_name_tot(iiko)=[];
    handles.tag_list_label.Value(iiko) = [];
    for j = 1:length(lab_name_tot)
        lab_name = lab_name_tot{j};
        if isequal(lab_name(1:6),'<HTML>');
            ii1 = find(lab_name=='>');
            ii2 = find(lab_name=='[');
            lab_name = lab_name(ii1(2)+1:ii2(1)-2);
        end
        
        xyz = handles.Labels_anatomy.(lab_name).xyz;
        hold(handles.tag_mesh,'on');
        aa(j) = plot3(xyz(:,1),xyz(:,2),xyz(:,3),'diamondk','markerfacecolor',cmap(j,:));
        map(handles.Labels_anatomy.(lab_name).index) = j;
    end
    
    h = findobj(handles.tag_mesh,'type','patch','tag','Mesh_Ventricles');
    h.FaceVertexCData = map;
    set(h,'facecolor','interp')
    cmap = lines(25);
    cmap(1,:) = [1 1 1]*.8;
    %     cmap=cmap(randperm(size(cmap,1)),:);
    colormap(handles.tag_mesh,cmap);
    caxis([0 22])
    
    
else
    cc = findobj(handles.tag_mesh,'marker','diamond');
    delete(cc);
    h = findobj(handles.tag_mesh,'type','patch','tag','Mesh_Ventricles');
    set(h,'edgecolor','k','marker','o','markeredgecolor','k','markerfacecolor','flat','facecolor',[1 1 1]*.7,'markerfacecolor','w','markersize',6);
    
end


% --- Executes on button press in tag_modify.
function tag_modify_Callback(hObject, eventdata, handles)
% hObject    handle to tag_modify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cc = findobj(handles.tag_mesh,'marker','diamond');
delete(cc);

if length(handles.tag_list_label.Value)>1
    handles.tag_list_label.Value = handles.tag_list_label.Value(1);
end
lab_name = handles.tag_list_label.String{handles.tag_list_label.Value};
if isequal(lab_name(1:6),'<HTML>');
    ii1 = find(lab_name=='>');
    ii2 = find(lab_name=='[');
    lab_name = lab_name(ii1(2)+1:ii2(1)-2);
end

xyz = handles.Labels_anatomy.(lab_name).xyz;
hold(handles.tag_mesh,'on');
for i=1:size(xyz,1)
    aa(i) = plot3(xyz(i,1),xyz(i,2),xyz(i,3),'or','markerfacecolor','g','markersize',8);
    aa(i).Tag = 'segment';
end


% --- Executes on button press in tag_delete.
function tag_delete_Callback(hObject, eventdata, handles)
% hObject    handle to tag_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lab_name_tot = handles.tag_list_label.String(handles.tag_list_label.Value);
for j = 1:length(lab_name_tot)
    lab_name = lab_name_tot{j};
    if isequal(lab_name(1:6),'<HTML>');
        ii1 = find(lab_name=='>');
        ii2 = find(lab_name=='[');
        lab_name = lab_name(ii1(2)+1:ii2(1)-2);
    end
    if ~isfield(handles.Labels_anatomy,lab_name);
        continue
    end
    A = questdlg(['Do you want to delete ',lab_name,' ?']);
    if isequal(A,'Yes')
        handles.Labels_anatomy = rmfield(handles.Labels_anatomy,lab_name);
        iistr = handles.tag_list_label.Value(j);
        str = handles.tag_list_label.String;
        str{iistr} = lab_name;
        handles.tag_list_label.String = str;
    end
end
% Update handles structure
cc = findobj(handles.tag_mesh,'marker','diamond');
delete(cc);
% h = findobj(handles.tag_mesh,'type','patch','tag','Mesh_Ventricles');
% set(h,'edgecolor','k','marker','o','markeredgecolor','k','markerfacecolor','flat','facecolor',[1 1 1]*.7,'markerfacecolor','w','markersize',6);

handles.tag_list_label.Value = 1:length(handles.tag_list_label.String);
tag_show_roi_Callback(hObject, [], handles);

% updating number of nodes
ii_existing = [];
vv = fieldnames(handles.Labels_anatomy);
for j = 1:length(vv)
    ii_existing = [ii_existing;handles.Labels_anatomy.(vv{j}).index];
end
nn = setdiff([1:size(handles.handles_ori.matfile.geo.Nodes_ventricles,1)],ii_existing);      
handles.tag_n_nodes.String = num2str(length(nn));

guidata(hObject, handles);


% --- Executes on button press in tag_overwrite.
function tag_overwrite_Callback(hObject, eventdata, handles)
% hObject    handle to tag_overwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_overwrite


% --- Executes on button press in tag_done.
function tag_done_Callback(hObject, eventdata, handles)
% hObject    handle to tag_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.handles_ori.matfile.geo.Labels_anatomy = handles.Labels_anatomy;
handles.tag_done.BackgroundColor = [1 0 0];

guidata(handles.handles_ori.figure1,handles.handles_ori);
A = questdlg('Did you finish the segmentation?');
% pause(0.25)
handles.tag_done.BackgroundColor = [1 1 1]*.94;
if isequal(A,'Yes');
    close(handles.figure1);
end


% --- Executes on button press in tag_spare.
function tag_spare_Callback(hObject, eventdata, handles)
% hObject    handle to tag_spare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vv = fieldnames(handles.Labels_anatomy);
ii = [];
for i = 1:length(vv)
    ii = [ii;handles.Labels_anatomy.(vv{i}).index];
    
end

xyz0 = handles.handles_ori.matfile.geo.Nodes_ventricles;

ii_spare = setdiff([1:size(xyz0,1)],ii);

cc1 = findobj(handles.tag_mesh,'markerfacecolor','g');
cc2 = findobj(handles.tag_mesh,'markerfacecolor','k');
cc3 = findobj(handles.tag_mesh,'markerfacecolor','r');
delete([cc1 cc2 cc3]);


aa = gobjects(length(ii_spare),1);
hold(handles.tag_mesh,'on');
for i=1:length(ii_spare)
    aa(i) = plot3(xyz0(ii_spare(i),1),xyz0(ii_spare(i),2),xyz0(ii_spare(i),3),'or','markerfacecolor','g','markersize',8);
    aa(i).Tag = 'segment';
end

handles.tag_show_roi.Value = 1;
handles.tag_list_label.Value = 1:length(handles.tag_list_label.String);
tag_show_roi_Callback(hObject, [], handles);



function tag_n_nodes_Callback(hObject, eventdata, handles)
% hObject    handle to tag_n_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tag_n_nodes as text
%        str2double(get(hObject,'String')) returns contents of tag_n_nodes as a double


% --- Executes during object creation, after setting all properties.
function tag_n_nodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tag_n_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
