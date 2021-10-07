function varargout = GUI_measure_distance(varargin)
% GUI_MEASURE_DISTANCE MATLAB code for GUI_measure_distance.fig
%      GUI_MEASURE_DISTANCE, by itself, creates a new GUI_MEASURE_DISTANCE or raises the existing
%      singleton*.
%
%      H = GUI_MEASURE_DISTANCE returns the handle to a new GUI_MEASURE_DISTANCE or the handle to
%      the existing singleton*.
%
%      GUI_MEASURE_DISTANCE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MEASURE_DISTANCE.M with the given input arguments.
%
%      GUI_MEASURE_DISTANCE('Property','Value',...) creates a new GUI_MEASURE_DISTANCE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_measure_distance_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_measure_distance_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_measure_distance

% Last Modified by GUIDE v2.5 28-Feb-2018 17:07:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_measure_distance_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_measure_distance_OutputFcn, ...
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


% --- Executes just before GUI_measure_distance is made visible.
function GUI_measure_distance_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_measure_distance (see VARARGIN)

% Choose default command line output for GUI_measure_distance
handles.output = hObject;

handles.figure_map = varargin{1};
dcm_obj = datacursormode(handles.figure_map);

%%
set(handles.tag_table_dist,'CellSelectionCallback',{@myfun_CellSelectionCallback_meas_dist})
handles.tag_table_dist.Data = cell(10,3);
handles.Points = cell(size(handles.tag_table_dist.Data));
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_measure_distance wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_measure_distance_OutputFcn(hObject, eventdata, handles)
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

h_main = guidata(handles.figure_map);


if handles.tag_add_from_file.Value==1
    [FILENAME, PATHNAME] = uigetfile('*Positions.txt','MultiSelect', 'off');
    
    loadname_pos = [PATHNAME,FILENAME];
    
    ii1 = findstr(FILENAME,'_P');ii1 = ii1(1);
    ii2 = findstr(FILENAME,'_');ii2(ii2<=ii1)=[];ii2 = ii2(1);
    pnumber = FILENAME(ii1+1:ii2-1);
    clear ii1 ii2
    
    if ~isempty(strfind(loadname_pos,'NAVISTAR_CONNECTOR_Eleclectrode_Positions.txt'))
        ll = {'M1','M2','M3','M4'};
        chan =  listdlg('liststring',ll,'PromptString',['Select ELECTRODES for P=',pnumber],'SelectionMode','single');
    elseif ~isempty(strfind(loadname_pos,'_20_POLE_A_CONNECTOR_Eleclectrode_Positions.txt'))|~isempty(strfind(loadname_pos,'_20_POLE_B_CONNECTOR_Eleclectrode_Positions.txt'))
        ll = {'Penta-1','Penta-2','Penta-3','Penta-4','Penta-5','Penta-6','Penta-7','Penta-8','Penta-9','Penta-10',...
            'Penta-11','Penta-12','Penta-13','Penta-14','Penta-15','Penta-16','Penta-17','Penta-18','Penta-19','Penta-20'};
        chan =  listdlg('liststring',ll,'PromptString',['Select ELECTRODES for P=',pnumber],'SelectionMode','single');
    else
        error('Point should be either a mapping or a penta-array position point');
    end
    
    Pos = LoadElPositionFromCarto(loadname_pos);
    
    hold(h_main.tag_mesh,'on');
    plot3(h_main.tag_mesh,Pos.xyz_OnAnnot(chan,1),Pos.xyz_OnAnnot(chan,2),Pos.xyz_OnAnnot(chan,3),'x','markersize',20,'linewidth',2,'color','k');
    
    i1 = handles.ii_select(1);
    i2 = handles.ii_select(2);
    
    handles.tag_table_dist.Data{i1,i2} = ['Ext.',num2str(chan)];
    handles.Points{i1,i2} = Pos.xyz_OnAnnot(chan,:);
else
    
    dcm_obj = datacursormode(handles.figure_map);
    info_struct = getCursorInfo(dcm_obj);
    
    [~,im] = min(sum(abs(h_main.matfile.geo.xyz  - repmat(info_struct.Position,[size(h_main.matfile.geo.xyz,1),1])),2));
    
    plot3(h_main.tag_mesh,h_main.matfile.geo.xyz(im,1),h_main.matfile.geo.xyz(im,2),h_main.matfile.geo.xyz(im,3),'*k','markersize',10);
    
    %     T = handles.tag_table_dist.Data;%
    i1 = handles.ii_select(1);
    i2 = handles.ii_select(2);
    %     T{i1,i2} = im;
    handles.tag_table_dist.Data{i1,i2} = im;
    handles.Points{i1,i2} = h_main.matfile.geo.xyz(im,:);
    
end
guidata(hObject,handles);

if ~isempty(handles.Points{i1,1})&~isempty(handles.Points{i1,2})
    p1 = handles.Points{i1,1};
    p2 = handles.Points{i1,2};
    D = sqrt(sum((p1-p2).^2,2));
    handles.tag_table_dist.Data{i1,3} = D;
    plot3(h_main.tag_mesh,[p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'-r');
    
end

% --- Executes on button press in tag_add_from_file.
function tag_add_from_file_Callback(hObject, eventdata, handles)
% hObject    handle to tag_add_from_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tag_add_from_file


% --- Executes on button press in tag_clear_map.
function tag_clear_map_Callback(hObject, eventdata, handles)
% hObject    handle to tag_clear_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

datacursormode off;
h_main = guidata(handles.figure_map);
mm = findobj(h_main.tag_mesh,'marker','x');
delete(mm);
mm = findobj(h_main.tag_mesh,'marker','*');
delete(mm);
mm = findobj(h_main.tag_mesh,'linestyle','-','color','r');
delete(mm);


% --- Executes on button press in tag_export_Table.
function tag_export_Table_Callback(hObject, eventdata, handles)
% hObject    handle to tag_export_Table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FILENAME, PATHNAME] = uiputfile('*.xls');

T1 = handles.tag_table_dist.Data;
if size(T1,2)==3
    T = [{'P1 (ind)','P2 (ind)','Dist (mm)'};T1];
elseif size(T1,2)==4
    T = [{'P1 (ind)','P2 (ind)','Dist (mm)','Comments'};T1];
else
    disp('Complete table before saving')
    return
end

T2 = handles.Points;
T2b = cell(size(T2,1),6);
if size(T1,2)>=2
    for j = 1:size(T2,1)
        if ~isempty(T2{j,1})
            T2b(j,1:3) = num2cell(T2{j,1});
        end
        if ~isempty(T2{j,2})
            T2b(j,4:6) = num2cell(T2{j,2});
        end
    end
else
    disp('Complete table before saving')
    return
end


T2b = [{'x','y','z','x','y','z'};T2b];

disp(['Writing ',PATHNAME,FILENAME]);
xlswrite([PATHNAME,FILENAME],T2b,'xyz','A1');
xlswrite([PATHNAME,FILENAME],T,'Dist','A1');
